#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/foreach_nonzero_with_fill_blocked.hpp>
#include <ctldl/symbolic/lower_triangle_blocked.hpp>

#include <cstddef>
#include <vector>

namespace ctldl {

constexpr auto getFilledInSparsityBlocked3x3(
    const SparsityView sparsity11_in, const SparsityView sparsity21_in,
    const SparsityView sparsity22_in, const SparsityView sparsity31_in,
    const SparsityView sparsity32_in, const SparsityView sparsity33_in,
    const PermutationView permutation1, const PermutationView permutation2,
    const PermutationView permutation3) {
  pre(isSquare(sparsity11_in));
  pre(isSquare(sparsity22_in));
  pre(sparsity21_in.numCols() == sparsity11_in.numCols());
  pre(sparsity22_in.numRows() == sparsity21_in.numRows());

  const auto sparsity11 = SparsityDynamicCSR(
      getSparsityDynamicLowerTriangle(sparsity11_in, permutation1));
  const auto sparsity21 = SparsityDynamicCSR(
      getSparsityDynamicPermuted(sparsity21_in, permutation2, permutation1));
  const auto sparsity22 = SparsityDynamicCSR(
      getSparsityDynamicLowerTriangle(sparsity22_in, permutation2));
  const auto sparsity31 = SparsityDynamicCSR(
      getSparsityDynamicPermuted(sparsity31_in, permutation3, permutation1));
  const auto sparsity32 = SparsityDynamicCSR(
      getSparsityDynamicPermuted(sparsity32_in, permutation3, permutation2));
  const auto sparsity33 = SparsityDynamicCSR(
      getSparsityDynamicLowerTriangle(sparsity33_in, permutation3));

  std::vector<Entry> entries11;
  std::vector<Entry> entries21;
  std::vector<Entry> entries22;
  std::vector<Entry> entries31;
  std::vector<Entry> entries32;
  std::vector<Entry> entries33;
  const auto add_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (j < sparsity11.numCols()) {
      if (i < sparsity11.numRows()) {
        entries11.push_back(Entry{i, j});
      } else if (i < sparsity11.numRows() + sparsity22.numRows()) {
        entries21.push_back(Entry{i - sparsity11.numRows(), j});
      } else {
        entries31.push_back(
            Entry{i - sparsity11.numRows() - sparsity22.numRows(), j});
      }
    } else if (j < sparsity11.numCols() + sparsity22.numCols()) {
      if (i < sparsity11.numRows() + sparsity22.numRows()) {
        entries22.push_back(
            Entry{i - sparsity11.numRows(), j - sparsity11.numCols()});
      } else {
        entries32.push_back(
            Entry{i - sparsity11.numRows() - sparsity22.numRows(),
                  j - sparsity11.numCols()});
      }
    } else {
      entries33.push_back(
          Entry{i - sparsity11.numRows() - sparsity22.numRows(),
                j - sparsity11.numCols() - sparsity22.numCols()});
    }
  };

  foreachNonZeroWithFillBlocked3x3(sparsity11, sparsity21, sparsity22,
                                   sparsity31, sparsity32, sparsity33,
                                   add_nonzero);
  // sorting is not needed for correctness, but helps performance
  sortEntriesRowMajorSorted(entries11);
  sortEntriesRowMajorSorted(entries21);
  sortEntriesRowMajorSorted(entries22);
  sortEntriesRowMajorSorted(entries31);
  sortEntriesRowMajorSorted(entries32);
  sortEntriesRowMajorSorted(entries33);

  return LowerTriangleBlocked3x3{
      SparsityDynamic(sparsity11.numRows(), sparsity11.numCols(), entries11),
      SparsityDynamic(sparsity21.numRows(), sparsity21.numCols(), entries21),
      SparsityDynamic(sparsity22.numRows(), sparsity22.numCols(), entries22),
      SparsityDynamic(sparsity31.numRows(), sparsity31.numCols(), entries31),
      SparsityDynamic(sparsity32.numRows(), sparsity32.numCols(), entries32),
      SparsityDynamic(sparsity33.numRows(), sparsity33.numCols(), entries33)};
}

constexpr auto getFilledInSparsityBlocked(const SparsityView sparsity11,
                                          const SparsityView sparsity21,
                                          const SparsityView sparsity22,
                                          const PermutationView permutation1,
                                          const PermutationView permutation2) {
  const auto sparsity31_dummy =
      makeEmptySparsityDynamic(0, sparsity11.numCols());
  const auto sparsity32_dummy =
      makeEmptySparsityDynamic(0, sparsity22.numCols());
  const auto sparsity33_dummy = makeEmptySparsityDynamic(0, 0);
  const auto permutation3_dummy = PermutationDynamic(0);
  const auto sparsity = getFilledInSparsityBlocked3x3(
      sparsity11, sparsity21, sparsity22, sparsity31_dummy, sparsity32_dummy,
      sparsity33_dummy, permutation1, permutation2, permutation3_dummy);
  return LowerTriangleBlocked{sparsity.block11, sparsity.block21,
                              sparsity.block22};
}

constexpr auto getFilledInSparsityBlocked(const SparsityView sparsity11,
                                          const SparsityView sparsity21,
                                          const SparsityView sparsity22) {
  return getFilledInSparsityBlocked(sparsity11, sparsity21, sparsity22,
                                    PermutationDynamic(sparsity11.numRows()),
                                    PermutationDynamic(sparsity22.numRows()));
}

}  // namespace ctldl
