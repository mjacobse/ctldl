#pragma once

#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permuted_entry_lower_triangle.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>
#include <ctldl/utility/contracts.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

constexpr auto getNnzLowerTriangle(const SparsityView sparsity) {
  pre(isSquare(sparsity));
  std::size_t nnz = 0;
  for (const auto entry : sparsity.entries()) {
    nnz += (entry.row_index > entry.col_index);
  }
  return nnz;
}

constexpr auto getEntriesLowerTriangle(const SparsityView sparsity,
                                       const PermutationView permutation) {
  pre(isSquare(sparsity));

  const auto inverse_permutation = invertPermutation(permutation);
  const auto is_lower_triangle = [](const Entry entry) {
    return entry.row_index > entry.col_index;
  };
  const auto permute_entry = [&inverse_permutation](const Entry entry) {
    return permutedEntryLowerTriangle(entry, inverse_permutation);
  };
  auto entries_permuted_lower_triangle = sparsity.entries() |
                                         std::views::filter(is_lower_triangle) |
                                         std::views::transform(permute_entry);
  return std::vector<Entry>(entries_permuted_lower_triangle.begin(),
                            entries_permuted_lower_triangle.end());
}

constexpr auto getSparsityDynamicLowerTriangle(
    const SparsityView sparsity, const PermutationView permutation) {
  pre(isSquare(sparsity));
  return SparsityDynamic(sparsity.numRows(), sparsity.numCols(),
                         getEntriesLowerTriangle(sparsity, permutation));
}

template <SparsityStatic sparsity>
constexpr auto getSparsityStaticLowerTriangle(
    const PermutationView permutation) {
  constexpr auto nnz = getNnzLowerTriangle(sparsity);
  std::array<Entry, nnz> entries;
  fixInitIfZeroLengthArray(entries);
  std::ranges::copy(getEntriesLowerTriangle(sparsity, permutation),
                    entries.begin());
  return SparsityStatic<nnz, sparsity.numRows(), sparsity.numCols()>(entries);
}

}  // namespace ctldl
