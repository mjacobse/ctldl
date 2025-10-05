#pragma once

#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>
#include <optional>
#include <vector>

namespace ctldl {

constexpr bool isSparsitySubset(
    const SparsityView sparsity_lhs_in, const SparsityView sparsity_rhs_in,
    const std::optional<PermutationView> permutation_row = std::nullopt,
    const std::optional<PermutationView> permutation_col = std::nullopt) {
  pre(sparsity_lhs_in.numRows() == sparsity_rhs_in.numRows());
  pre(sparsity_lhs_in.numCols() == sparsity_rhs_in.numCols());

  const auto sparsity_lhs = SparsityDynamicCSR(getSparsityDynamicPermuted(
      sparsity_lhs_in,
      permutation_row.value_or(PermutationDynamic(sparsity_lhs_in.numRows())),
      permutation_col.value_or(PermutationDynamic(sparsity_lhs_in.numCols()))));
  const auto sparsity_rhs = SparsityDynamicCSR(
      SparsityDynamic(sparsity_rhs_in.numRows(), sparsity_rhs_in.numCols(),
                      sparsity_rhs_in.entries()));

  const auto num_rows = std::size_t{sparsity_lhs.numRows()};
  const auto num_cols = std::size_t{sparsity_lhs.numCols()};

  std::vector<bool> is_nonzero_rhs(num_cols, false);
  for (std::size_t row_index = 0; row_index < num_rows; ++row_index) {
    for (const auto entry : sparsity_rhs.rowView(row_index)) {
      is_nonzero_rhs[entry.col_index] = true;
    }

    for (const auto entry : sparsity_lhs.rowView(row_index)) {
      if (!is_nonzero_rhs[entry.col_index]) {
        return false;
      }
    }

    // undo marked entries to start off clean in next row
    for (const auto entry : sparsity_rhs.rowView(row_index)) {
      is_nonzero_rhs[entry.col_index] = false;
    }
  }

  return true;
}

constexpr bool isSparsitySubsetLowerTriangle(
    const SparsityView sparsity_lhs, const SparsityView sparsity_rhs,
    const PermutationView permutation) {
  return isSparsitySubset(
      getSparsityDynamicLowerTriangle(sparsity_lhs, permutation), sparsity_rhs);
}

}  // namespace ctldl
