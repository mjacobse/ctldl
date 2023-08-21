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

template <class Sparsity11, class Sparsity21, class Sparsity22>
constexpr auto getFilledInNumNonZerosBlocked(const Sparsity11& sparsity11,
                                             const Sparsity21& sparsity21,
                                             const Sparsity22& sparsity22) {
  static_assert(isSquare<Sparsity11>());
  static_assert(isSquare<Sparsity22>());
  static_assert(Sparsity21::num_cols == Sparsity11::num_cols);
  static_assert(Sparsity22::num_rows == Sparsity21::num_rows);

  LowerTriangleBlocked nnz{std::size_t{0}, std::size_t{0}, std::size_t{0}};
  const auto count_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (j < sparsity11.num_cols) {
      if (i < sparsity11.num_rows) {
        nnz.block11 += 1;
      } else {
        nnz.block21 += 1;
      }
    } else {
      nnz.block22 += 1;
    }
  };
  foreachNonZeroWithFillBlocked(sparsity11, sparsity21, sparsity22,
                                count_nonzero);
  return nnz;
}

template <auto sparsity11_in, auto sparsity21_in, auto sparsity22_in,
          auto permutation1 = Permutation<sparsity11_in.num_rows>(),
          auto permutation2 = Permutation<sparsity22_in.num_rows>()>
constexpr auto getFilledInSparsityBlocked() {
  static_assert(isSquare(sparsity11_in));
  static_assert(isSquare(sparsity22_in));
  static_assert(sparsity21_in.num_cols == sparsity11_in.num_cols);
  static_assert(sparsity22_in.num_rows == sparsity21_in.num_rows);

  constexpr auto sparsity11 =
      SparsityCSR(getSparsityLowerTriangle<sparsity11_in>(permutation1));
  constexpr auto sparsity21 = SparsityCSR(
      getSparsityPermuted(sparsity21_in, permutation2, permutation1));
  constexpr auto sparsity22 =
      SparsityCSR(getSparsityLowerTriangle<sparsity22_in>(permutation2));

  constexpr auto nnz =
      getFilledInNumNonZerosBlocked(sparsity11, sparsity21, sparsity22);

  std::array<Entry, nnz.block11> entries11;
  std::array<Entry, nnz.block21> entries21;
  std::array<Entry, nnz.block22> entries22;
  LowerTriangleBlocked entry_index{std::size_t{0}, std::size_t{0},
                                   std::size_t{0}};
  const auto add_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (j < sparsity11.num_cols) {
      if (i < sparsity11.num_rows) {
        entries11[entry_index.block11] = Entry{i, j};
        entry_index.block11 += 1;
      } else {
        entries21[entry_index.block21] = Entry{i - sparsity11.num_rows, j};
        entry_index.block21 += 1;
      }
    } else {
      entries22[entry_index.block22] =
          Entry{i - sparsity11.num_rows, j - sparsity11.num_cols};
      entry_index.block22 += 1;
    }
  };

  foreachNonZeroWithFillBlocked(sparsity11, sparsity21, sparsity22,
                                add_nonzero);
  // sorting is not needed for correctness, but helps performance
  sortEntriesRowMajorSorted(entries11);
  sortEntriesRowMajorSorted(entries21);
  sortEntriesRowMajorSorted(entries22);

  return LowerTriangleBlocked{
      makeSparsity<sparsity11.num_rows, sparsity11.num_cols>(entries11),
      makeSparsity<sparsity21.num_rows, sparsity21.num_cols>(entries21),
      makeSparsity<sparsity22.num_rows, sparsity22.num_cols>(entries22)};
}

}  // namespace ctldl
