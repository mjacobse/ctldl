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

#include <array>
#include <cstddef>

namespace ctldl {

template <class SparsityDiag, class SparsitySubdiag, class SparsityOuter>
constexpr auto getNumNonZerosWithFillRepeatingArrowhead(
    const SparsityDiag& sparsity_A, const SparsitySubdiag& sparsity_B,
    const SparsityOuter& sparsity_C) {
  static_assert(isSquare<SparsityDiag>());
  static_assert(isSquare<SparsitySubdiag>());
  static_assert(SparsityDiag::numRows() == SparsitySubdiag::numRows());
  constexpr auto dim = std::size_t{SparsityDiag::numRows()};
  static_assert(SparsityOuter::numCols() == SparsityDiag::numCols());

  RepeatingBlockTridiagonalArrowhead nnz{std::size_t{0}, std::size_t{0},
                                         std::size_t{0}};
  const auto count_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (i >= 2 * dim) {
      nnz.outer += 1;
    } else {
      if (j >= dim) {
        nnz.diag += 1;
      } else {
        nnz.subdiag += 1;
      }
    }
  };
  foreachNonZeroWithFillRepeatingArrowhead(sparsity_A, sparsity_B, sparsity_C,
                                           count_nonzero);
  return nnz;
}

template <auto sparsity_in_A, auto sparsity_in_B, auto sparsity_in_C,
          auto permutation, auto permutation_C>
constexpr auto getFilledInSparsityRepeatingArrowhead() {
  static_assert(isSquare(sparsity_in_A));
  static_assert(isSquare(sparsity_in_B));
  static_assert(sparsity_in_B.numRows() == sparsity_in_A.numRows());
  static_assert(sparsity_in_C.numCols() == sparsity_in_A.numCols());
  constexpr auto dim = std::size_t{sparsity_in_A.numRows()};
  constexpr auto dim_outer = std::size_t{sparsity_in_C.numRows()};

  constexpr auto sparsity_A =
      SparsityCSR(getSparsityLowerTriangle<sparsity_in_A>(permutation));
  constexpr auto sparsity_B =
      SparsityCSR(getSparsityPermuted(sparsity_in_B, permutation, permutation));
  constexpr auto sparsity_C = SparsityCSR(
      getSparsityPermuted(sparsity_in_C, permutation_C, permutation));

  constexpr auto nnz = getNumNonZerosWithFillRepeatingArrowhead(
      sparsity_A, sparsity_B, sparsity_C);

  RepeatingBlockTridiagonalArrowhead entries{std::array<Entry, nnz.diag>{},
                                             std::array<Entry, nnz.subdiag>{},
                                             std::array<Entry, nnz.outer>{}};
  RepeatingBlockTridiagonalArrowhead entry_index{std::size_t{0}, std::size_t{0},
                                                 std::size_t{0}};
  const auto add_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (i >= 2 * dim) {
      entries.outer[entry_index.outer] = Entry{i - 2 * dim, j % dim};
      entry_index.outer += 1;
    } else {
      if (j >= dim) {
        entries.diag[entry_index.diag] = Entry{i - dim, j - dim};
        entry_index.diag += 1;
      } else {
        entries.subdiag[entry_index.subdiag] = Entry{i - dim, j};
        entry_index.subdiag += 1;
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
      makeSparsity<dim, dim>(entries.diag),
      makeSparsity<dim, dim>(entries.subdiag),
      makeSparsity<dim_outer, dim>(entries.outer)};
};

}  // namespace ctldl
