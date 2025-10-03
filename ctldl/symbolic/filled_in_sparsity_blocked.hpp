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

#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity11, class Sparsity21, class Sparsity22,
          class Sparsity31, class Sparsity32, class Sparsity33>
constexpr auto getFilledInNumNonZerosBlocked3x3(const Sparsity11& sparsity11,
                                                const Sparsity21& sparsity21,
                                                const Sparsity22& sparsity22,
                                                const Sparsity31& sparsity31,
                                                const Sparsity32& sparsity32,
                                                const Sparsity33& sparsity33) {
  LowerTriangleBlocked3x3 nnz{std::size_t{0}, std::size_t{0}, std::size_t{0},
                              std::size_t{0}, std::size_t{0}, std::size_t{0}};
  const auto count_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (j < sparsity11.numCols()) {
      if (i < sparsity11.numRows()) {
        nnz.block11 += 1;
      } else if (i < sparsity11.numRows() + sparsity22.numRows()) {
        nnz.block21 += 1;
      } else {
        nnz.block31 += 1;
      }
    } else if (j < sparsity11.numCols() + sparsity22.numCols()) {
      if (i < sparsity11.numRows() + sparsity22.numRows()) {
        nnz.block22 += 1;
      } else {
        nnz.block32 += 1;
      }
    } else {
      nnz.block33 += 1;
    }
  };
  foreachNonZeroWithFillBlocked3x3(sparsity11, sparsity21, sparsity22,
                                   sparsity31, sparsity32, sparsity33,
                                   count_nonzero);
  return nnz;
}

template <class Sparsity11, class Sparsity21, class Sparsity22>
constexpr auto getFilledInNumNonZerosBlocked(const Sparsity11& sparsity11,
                                             const Sparsity21& sparsity21,
                                             const Sparsity22& sparsity22) {
  constexpr auto sparsity31_dummy =
      makeEmptySparsityCSR<0, Sparsity11::numCols()>();
  constexpr auto sparsity32_dummy =
      makeEmptySparsityCSR<0, Sparsity22::numCols()>();
  constexpr auto sparsity33_dummy = makeEmptySparsityCSR<0, 0>();
  const auto nnz = getFilledInNumNonZerosBlocked3x3(
      sparsity11, sparsity21, sparsity22, sparsity31_dummy, sparsity32_dummy,
      sparsity33_dummy);
  return LowerTriangleBlocked{nnz.block11, nnz.block21, nnz.block22};
}

template <auto sparsity11_in, auto sparsity21_in, auto sparsity22_in,
          auto sparsity31_in, auto sparsity32_in, auto sparsity33_in,
          auto permutation1 = Permutation<sparsity11_in.numRows()>(),
          auto permutation2 = Permutation<sparsity22_in.numRows()>(),
          auto permutation3 = Permutation<sparsity33_in.numRows()>()>
constexpr auto getFilledInSparsityBlocked3x3() {
  static_assert(isSquare(sparsity11_in));
  static_assert(isSquare(sparsity22_in));
  static_assert(sparsity21_in.numCols() == sparsity11_in.numCols());
  static_assert(sparsity22_in.numRows() == sparsity21_in.numRows());

  constexpr auto sparsity11 =
      SparsityCSR(getSparsityLowerTriangle<sparsity11_in>(permutation1));
  constexpr auto sparsity21 = SparsityCSR(
      getSparsityPermuted(sparsity21_in, permutation2, permutation1));
  constexpr auto sparsity22 =
      SparsityCSR(getSparsityLowerTriangle<sparsity22_in>(permutation2));
  constexpr auto sparsity31 = SparsityCSR(
      getSparsityPermuted(sparsity31_in, permutation3, permutation1));
  constexpr auto sparsity32 = SparsityCSR(
      getSparsityPermuted(sparsity32_in, permutation3, permutation2));
  constexpr auto sparsity33 =
      SparsityCSR(getSparsityLowerTriangle<sparsity33_in>(permutation3));

  constexpr auto nnz = getFilledInNumNonZerosBlocked3x3(
      sparsity11, sparsity21, sparsity22, sparsity31, sparsity32, sparsity33);

  std::array<Entry, nnz.block11> entries11;
  std::array<Entry, nnz.block21> entries21;
  std::array<Entry, nnz.block22> entries22;
  std::array<Entry, nnz.block31> entries31;
  std::array<Entry, nnz.block32> entries32;
  std::array<Entry, nnz.block33> entries33;
  LowerTriangleBlocked3x3 entry_index{std::size_t{0}, std::size_t{0},
                                      std::size_t{0}, std::size_t{0},
                                      std::size_t{0}, std::size_t{0}};
  const auto add_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (j < sparsity11.numCols()) {
      if (i < sparsity11.numRows()) {
        entries11[entry_index.block11] = Entry{i, j};
        entry_index.block11 += 1;
      } else if (i < sparsity11.numRows() + sparsity22.numRows()) {
        entries21[entry_index.block21] = Entry{i - sparsity11.numRows(), j};
        entry_index.block21 += 1;
      } else {
        entries31[entry_index.block31] =
            Entry{i - sparsity11.numRows() - sparsity22.numRows(), j};
        entry_index.block31 += 1;
      }
    } else if (j < sparsity11.numCols() + sparsity22.numCols()) {
      if (i < sparsity11.numRows() + sparsity22.numRows()) {
        entries22[entry_index.block22] =
            Entry{i - sparsity11.numRows(), j - sparsity11.numCols()};
        entry_index.block22 += 1;
      } else {
        entries32[entry_index.block32] =
            Entry{i - sparsity11.numRows() - sparsity22.numRows(),
                  j - sparsity11.numCols()};
        entry_index.block32 += 1;
      }
    } else {
      entries33[entry_index.block33] =
          Entry{i - sparsity11.numRows() - sparsity22.numRows(),
                j - sparsity11.numCols() - sparsity22.numCols()};
      entry_index.block33 += 1;
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
      makeSparsity<sparsity11.numRows(), sparsity11.numCols()>(entries11),
      makeSparsity<sparsity21.numRows(), sparsity21.numCols()>(entries21),
      makeSparsity<sparsity22.numRows(), sparsity22.numCols()>(entries22),
      makeSparsity<sparsity31.numRows(), sparsity31.numCols()>(entries31),
      makeSparsity<sparsity32.numRows(), sparsity32.numCols()>(entries32),
      makeSparsity<sparsity33.numRows(), sparsity33.numCols()>(entries33)};
}

template <auto sparsity11, auto sparsity21, auto sparsity22,
          auto permutation1 = Permutation<sparsity11.numRows()>(),
          auto permutation2 = Permutation<sparsity22.numRows()>()>
constexpr auto getFilledInSparsityBlocked() {
  constexpr auto sparsity31_dummy = makeEmptySparsity<0, sparsity11.numCols()>();
  constexpr auto sparsity32_dummy = makeEmptySparsity<0, sparsity22.numCols()>();
  constexpr auto sparsity33_dummy = makeEmptySparsity<0, 0>();
  constexpr auto permutation3_dummy = Permutation<0>{};
  const auto sparsity = getFilledInSparsityBlocked3x3<
      sparsity11, sparsity21, sparsity22, sparsity31_dummy, sparsity32_dummy,
      sparsity33_dummy, permutation1, permutation2, permutation3_dummy>();
  return LowerTriangleBlocked{sparsity.block11, sparsity.block21,
                              sparsity.block22};
}

}  // namespace ctldl
