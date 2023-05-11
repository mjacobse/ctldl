#pragma once

#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity, class PermutationIn>
struct SparsityLowerTriangle {
  static_assert(Sparsity::num_rows == Sparsity::num_cols);
  static constexpr auto dim = std::size_t{Sparsity::num_rows};
  static constexpr auto num_rows = dim;
  static constexpr auto num_cols = dim;
  static constexpr auto nnz = [] {
    std::size_t nnz = 0;
    for (const auto entry : Sparsity::entries) {
      nnz += (entry.row_index > entry.col_index);
    }
    return nnz;
  }();

  static constexpr Permutation<dim> permutation{PermutationIn::permutation};
  static constexpr auto inverse_permutation = invertPermutation(permutation);

  static constexpr auto entries = [] {
    std::array<Entry, nnz> entries;
    std::size_t entry_index = 0;
    for (const auto entry : Sparsity::entries) {
      if (entry.row_index <= entry.col_index) {
        continue;
      }
      const auto new_row_index = std::max(inverse_permutation[entry.row_index],
                                          inverse_permutation[entry.col_index]);
      const auto new_col_index = std::min(inverse_permutation[entry.row_index],
                                          inverse_permutation[entry.col_index]);
      entries[entry_index] = Entry{new_row_index, new_col_index};
      entry_index += 1;
    }
    return entries;
  }();
};

}  // namespace ctldl
