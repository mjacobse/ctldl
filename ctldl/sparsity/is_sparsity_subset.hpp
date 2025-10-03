#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class SparsityLhs, class SparsityRhs>
constexpr bool isSparsitySubset(
    const SparsityLhs& sparsity_lhs_in, const SparsityRhs& sparsity_rhs_in,
    const PermutationView permutation_row =
        PermutationStatic<SparsityLhs::numRows()>{},
    const PermutationView permutation_col =
        PermutationStatic<SparsityLhs::numCols()>{}) {
  const auto sparsity_lhs = makeSparsityCSR(
      getSparsityPermuted(sparsity_lhs_in, permutation_row, permutation_col));
  const auto sparsity_rhs = makeSparsityCSR(sparsity_rhs_in);

  static_assert(sparsity_lhs.numRows() == sparsity_rhs.numRows());
  static_assert(sparsity_lhs.numCols() == sparsity_rhs.numCols());
  constexpr auto num_rows = std::size_t{sparsity_lhs.numRows()};
  constexpr auto num_cols = std::size_t{sparsity_lhs.numCols()};

  std::array<bool, num_cols> is_nonzero_rhs;
  fixInitIfZeroLengthArray(is_nonzero_rhs);
  is_nonzero_rhs.fill(false);
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

template <auto sparsity_lhs, class SparsityRhs>
constexpr bool isSparsitySubsetLowerTriangle(
    const SparsityRhs& sparsity_rhs, const PermutationView permutation) {
  return isSparsitySubset(getSparsityLowerTriangle<sparsity_lhs>(permutation),
                          sparsity_rhs);
}

}  // namespace ctldl
