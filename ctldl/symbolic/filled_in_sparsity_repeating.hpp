#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/foreach_nonzero_with_fill_repeated.hpp>
#include <ctldl/symbolic/repeating_block_tridiagonal.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class SparsityA, class SparsityB>
constexpr auto getNumNonZerosWithFillRepeating(const SparsityA& sparsity_A,
                                               const SparsityB& sparsity_B) {
  static_assert(isSquare<SparsityA>());
  static_assert(isSquare<SparsityB>());
  static_assert(SparsityA::num_rows == SparsityB::num_rows);
  constexpr auto dim = std::size_t{SparsityA::num_rows};

  auto nnz = RepeatingBlockTridiagonal{std::size_t{0}, std::size_t{0}};
  const auto count_nonzero = [&](const std::size_t /*i*/, const std::size_t j) {
    if (j >= dim) {
      nnz.diagonal += 1;
    } else {
      nnz.subdiagonal += 1;
    }
  };
  foreachNonZeroWithFillRepeated(sparsity_A, sparsity_B, count_nonzero);
  return nnz;
}

template <auto sparsity_in_A, auto sparsity_in_B, auto permutation>
constexpr auto getFilledInSparsityRepeating() {
  static_assert(isSquare(sparsity_in_A));
  static_assert(isSquare(sparsity_in_B));
  static_assert(sparsity_in_B.num_rows == sparsity_in_A.num_rows);
  constexpr auto dim = std::size_t{sparsity_in_A.num_rows};

  constexpr auto sparsity_A =
      SparsityCSR(getSparsityLowerTriangle<sparsity_in_A>(permutation));
  constexpr auto sparsity_B =
      SparsityCSR(getSparsityPermuted(sparsity_in_B, permutation, permutation));

  constexpr auto nnz = getNumNonZerosWithFillRepeating(sparsity_A, sparsity_B);

  std::array<Entry, nnz.diagonal> entries_A;
  std::array<Entry, nnz.subdiagonal> entries_B;
  RepeatingBlockTridiagonal entry_index{std::size_t{0}, std::size_t{0}};
  const auto add_nonzero = [&](const std::size_t i, const std::size_t j) {
    if (j >= dim) {
      entries_A[entry_index.diagonal] = Entry{i - dim, j - dim};
      entry_index.diagonal += 1;
    } else {
      entries_B[entry_index.subdiagonal] = Entry{i - dim, j};
      entry_index.subdiagonal += 1;
    }
  };
  foreachNonZeroWithFillRepeated(sparsity_A, sparsity_B, add_nonzero);

  // sorting is not needed for correctness, but helps performance
  sortEntriesRowMajorSorted(entries_A);
  sortEntriesRowMajorSorted(entries_B);

  return RepeatingBlockTridiagonal{makeSparsity<dim, dim>(entries_A),
                                   makeSparsity<dim, dim>(entries_B)};
};

}  // namespace ctldl
