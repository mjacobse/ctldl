#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class SparsityLhs, class SparsityRhs>
constexpr bool isSparsitySubset(
    const SparsityLhs& sparsity_lhs_in, const SparsityRhs& sparsity_rhs_in,
    const Permutation<SparsityLhs::num_rows>& permutation_row = {},
    const Permutation<SparsityLhs::num_cols>& permutation_col = {}) {
  const auto sparsity_lhs = makeSparsityCSR(
      getSparsityPermuted(sparsity_lhs_in, permutation_row, permutation_col));
  const auto sparsity_rhs = makeSparsityCSR(sparsity_rhs_in);

  static_assert(sparsity_lhs.num_rows == sparsity_rhs.num_rows);
  static_assert(sparsity_lhs.num_cols == sparsity_rhs.num_cols);
  constexpr auto num_rows = std::size_t{sparsity_lhs.num_rows};
  constexpr auto num_cols = std::size_t{sparsity_lhs.num_cols};

  std::array<bool, num_cols> is_nonzero_rhs;
  is_nonzero_rhs.fill(false);
  for (std::size_t row_index = 0; row_index < num_rows; ++row_index) {
    const auto row_begin_rhs = sparsity_rhs.row_begin_indices[row_index];
    const auto row_end_rhs = sparsity_rhs.row_begin_indices[row_index + 1];
    for (auto entry_index = row_begin_rhs; entry_index != row_end_rhs;
         ++entry_index) {
      const auto col_index = sparsity_rhs.entries[entry_index].col_index;
      is_nonzero_rhs[col_index] = true;
    }

    const auto row_begin_lhs = sparsity_lhs.row_begin_indices[row_index];
    const auto row_end_lhs = sparsity_lhs.row_begin_indices[row_index + 1];
    for (auto entry_index = row_begin_lhs; entry_index != row_end_lhs;
         ++entry_index) {
      const auto col_index = sparsity_lhs.entries[entry_index].col_index;
      if (!is_nonzero_rhs[col_index]) {
        return false;
      }
    }

    for (auto entry_index = row_begin_rhs; entry_index != row_end_rhs;
         ++entry_index) {
      const auto col_index = sparsity_rhs.entries[entry_index].col_index;
      is_nonzero_rhs[col_index] = false;
    }
  }

  return true;
}

template <auto sparsity_lhs, class SparsityRhs, class Permutation>
constexpr bool isSparsitySubsetLowerTriangle(const SparsityRhs& sparsity_rhs,
                                             const Permutation& permutation) {
  return isSparsitySubset(getSparsityLowerTriangle<sparsity_lhs>(permutation),
                          sparsity_rhs);
}

}  // namespace ctldl
