#pragma once

#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permuted_entry.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity, class PermutationRow, class PermutationCol>
struct SparsityPermuted {
  static constexpr auto num_rows = std::size_t{Sparsity::num_rows};
  static constexpr auto num_cols = std::size_t{Sparsity::num_cols};
  static constexpr auto nnz = std::size_t{std::size(Sparsity::entries)};

  static constexpr Permutation<num_rows> permutation_row{PermutationRow::permutation};
  static constexpr Permutation<num_cols> permutation_col{PermutationCol::permutation};
  static constexpr auto inverse_permutation_row = invertPermutation(permutation_row);
  static constexpr auto inverse_permutation_col = invertPermutation(permutation_col);

  static constexpr auto entries = [] {
    std::array<Entry, nnz> entries;
    std::size_t entry_index = 0;
    for (const auto entry : Sparsity::entries) {
      entries[entry_index] = permutedEntry(entry, inverse_permutation_row,
                                           inverse_permutation_col);
      entry_index += 1;
    }
    return entries;
  }();
};

}  // namespace ctldl
