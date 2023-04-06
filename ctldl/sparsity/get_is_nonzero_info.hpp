#pragma once

#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/is_nonzero_info.hpp>

#include <algorithm>

namespace ctldl {

template <class Sparsity>
constexpr auto getIsNonzeroInfo(
    const Permutation<Sparsity::num_rows>& permutation_row,
    const Permutation<Sparsity::num_cols>& permutation_col) {
  const auto inverse_permutation_row = invertPermutation(permutation_row);
  const auto inverse_permutation_col = invertPermutation(permutation_col);

  IsNonzeroInfo<Sparsity::num_rows, Sparsity::num_cols> is_nonzero;
  for (const auto entry : Sparsity::entries) {
    const auto new_row_index = inverse_permutation_row[entry.row_index];
    const auto new_col_index = inverse_permutation_col[entry.col_index];
    is_nonzero[new_row_index][new_col_index] = true;
  }
  return is_nonzero;
}

template <class Sparsity>
constexpr auto getIsNonzeroInfo() {
  return getIsNonzeroInfo<Sparsity>(PermutationIdentity::permutation,
                                    PermutationIdentity::permutation);
}

template <class Sparsity>
constexpr auto getIsNonzeroInfoLowerTriangle(
    const Permutation<Sparsity::num_rows>& permutation =
        PermutationIdentity::permutation) {
  static_assert(Sparsity::num_rows == Sparsity::num_cols);
  const auto inverse_permutation = invertPermutation(permutation);

  IsNonzeroInfo<Sparsity::num_rows, Sparsity::num_cols> is_nonzero;
  for (const auto entry : Sparsity::entries) {
    const auto new_row_index = std::max(inverse_permutation[entry.row_index],
                                        inverse_permutation[entry.col_index]);
    const auto new_col_index = std::min(inverse_permutation[entry.row_index],
                                        inverse_permutation[entry.col_index]);
    is_nonzero[new_row_index][new_col_index] = true;
  }
  return is_nonzero;
}

}  // namespace ctldl
