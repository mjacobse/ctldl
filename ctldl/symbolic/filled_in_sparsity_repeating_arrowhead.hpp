#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/foreach_nonzero_with_fill_repeating_arrowhead.hpp>
#include <ctldl/symbolic/repeating_block_tridiagonal_arrowhead.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>
#include <vector>

namespace ctldl {

constexpr auto getFilledInSparsityRepeatingArrowhead(
    const SparsityView sparsity_in_A, const SparsityView sparsity_in_B,
    const SparsityView sparsity_in_C, const PermutationView permutation,
    const PermutationView permutation_C) {
  pre(isSquare(sparsity_in_A));
  pre(isSquare(sparsity_in_B));
  pre(sparsity_in_B.numRows() == sparsity_in_A.numRows());
  pre(sparsity_in_C.numCols() == sparsity_in_A.numCols());
  const auto dim = std::size_t{sparsity_in_A.numRows()};
  const auto dim_outer = std::size_t{sparsity_in_C.numRows()};

  const auto sparsity_A = SparsityDynamicCSR(
      getSparsityDynamicLowerTriangle(sparsity_in_A, permutation));
  const auto sparsity_B = SparsityDynamicCSR(
      getSparsityDynamicPermuted(sparsity_in_B, permutation, permutation));
  const auto sparsity_C = SparsityDynamicCSR(
      getSparsityDynamicPermuted(sparsity_in_C, permutation_C, permutation));

  RepeatingBlockTridiagonalArrowhead entries{
      std::vector<Entry>{}, std::vector<Entry>{},
      std::vector<Entry>{}};
  const auto add_nonzero = [dim, &entries](const std::size_t i,
                                           const std::size_t j) {
    if (i >= 2 * dim) {
      entries.outer.push_back(Entry{i - 2 * dim, j % dim});
    } else {
      if (j >= dim) {
        entries.diag.push_back(Entry{i - dim, j - dim});
      } else {
        entries.subdiag.push_back(Entry{i - dim, j});
      }
    }
  };
  foreachNonZeroWithFillRepeatingArrowhead(sparsity_A, sparsity_B, sparsity_C,
                                           add_nonzero);

  // sorting is not needed for correctness, but helps performance
  sortEntriesRowMajorSorted(entries.diag);
  sortEntriesRowMajorSorted(entries.subdiag);
  sortEntriesRowMajorSorted(entries.outer);

  return RepeatingBlockTridiagonalArrowhead{
      SparsityDynamic(dim, dim, entries.diag),
      SparsityDynamic(dim, dim, entries.subdiag),
      SparsityDynamic(dim_outer, dim, entries.outer)};
};

}  // namespace ctldl
