#pragma once

#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permuted_entry.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <array>
#include <ranges>
#include <vector>

namespace ctldl {

constexpr auto getEntriesPermuted(const SparsityView sparsity,
                                  const PermutationView permutation_row,
                                  const PermutationView permutation_col) {
  const auto inverse_permutation_row = invertPermutation(permutation_row);
  const auto inverse_permutation_col = invertPermutation(permutation_col);

  const auto permute_entry = [&inverse_permutation_row,
                              &inverse_permutation_col](const Entry entry) {
    return permutedEntry(entry, inverse_permutation_row,
                         inverse_permutation_col);
  };
  auto entries_permuted =
      std::views::transform(sparsity.entries(), permute_entry);
  return std::vector<Entry>(entries_permuted.begin(), entries_permuted.end());
}

constexpr auto getSparsityDynamicPermuted(
    const SparsityView sparsity, const PermutationView permutation_row,
    const PermutationView permutation_col) {
  return SparsityDynamic(
      sparsity.numRows(), sparsity.numCols(),
      getEntriesPermuted(sparsity, permutation_row, permutation_col));
}

template <class SparsityIn>
constexpr auto getSparsityStaticPermuted(
    const SparsityIn& sparsity, const PermutationView permutation_row,
    const PermutationView permutation_col) {
  using std::size;
  constexpr auto nnz = std::size_t{size(decltype(sparsity.entries()){})};

  std::array<Entry, nnz> entries;
  std::ranges::copy(
      getEntriesPermuted(sparsity, permutation_row, permutation_col),
      entries.begin());
  return SparsityStatic<nnz, SparsityIn::numRows(), SparsityIn::numCols()>(
      entries);
}

}  // namespace ctldl
