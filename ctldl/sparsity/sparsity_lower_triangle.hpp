#pragma once

#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permuted_entry_lower_triangle.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity>
constexpr auto getNnzLowerTriangle(const Sparsity& sparsity) {
  static_assert(isSquare<Sparsity>());
  std::size_t nnz = 0;
  for (const auto entry : sparsity.entries()) {
    nnz += (entry.row_index > entry.col_index);
  }
  return nnz;
}

template <std::size_t nnz, class Sparsity, class PermutationIn>
constexpr auto getEntriesLowerTriangle(const Sparsity& sparsity,
                                       const PermutationIn& permutation) {
  static_assert(isSquare<Sparsity>());
  const auto inverse_permutation = invertPermutation(permutation);

  std::array<Entry, nnz> entries;
  fixInitIfZeroLengthArray(entries);
  std::size_t entry_index = 0;
  for (const auto entry : sparsity.entries()) {
    if (entry.row_index <= entry.col_index) {
      continue;
    }
    entries[entry_index] =
        permutedEntryLowerTriangle(entry, inverse_permutation);
    entry_index += 1;
  }
  return entries;
}

template <auto sparsity, class PermutationIn>
constexpr auto getSparsityLowerTriangle(const PermutationIn& permutation) {
  static_assert(isSquare(sparsity));
  constexpr auto dim = std::size_t{sparsity.numRows()};
  constexpr auto nnz = getNnzLowerTriangle(sparsity);
  const auto entries = getEntriesLowerTriangle<nnz>(sparsity, permutation);
  return Sparsity<nnz, dim, dim>(entries);
}

}  // namespace ctldl
