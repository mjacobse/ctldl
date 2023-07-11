#pragma once

#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permuted_entry.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class SparsityIn, class PermutationRow, class PermutationCol>
constexpr auto getSparsityPermuted(const SparsityIn& sparsity,
                                   const PermutationRow& permutation_row,
                                   const PermutationCol& permutation_col) {
  using std::size;
  constexpr auto nnz = std::size_t{size(decltype(sparsity.entries){})};

  const auto inverse_permutation_row = invertPermutation(permutation_row);
  const auto inverse_permutation_col = invertPermutation(permutation_col);

  std::array<Entry, nnz> entries;
  std::size_t entry_index = 0;
  for (const auto entry : sparsity.entries) {
    entries[entry_index] = permutedEntry(entry, inverse_permutation_row,
                                         inverse_permutation_col);
    entry_index += 1;
  }

  return Sparsity<nnz, SparsityIn::num_rows, SparsityIn::num_cols>(entries);
}

}  // namespace ctldl
