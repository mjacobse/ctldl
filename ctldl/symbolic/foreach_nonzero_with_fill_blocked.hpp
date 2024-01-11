#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/compute_elimination_tree_blocked.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity11, class Sparsity21, class Sparsity22,
          class Sparsity31, class Sparsity32, class Sparsity33,
          class UnaryFunction>
constexpr void foreachNonZeroWithFillBlocked3x3(const Sparsity11& sparsity11,
                                                const Sparsity21& sparsity21,
                                                const Sparsity22& sparsity22,
                                                const Sparsity31& sparsity31,
                                                const Sparsity32& sparsity32,
                                                const Sparsity33& sparsity33,
                                                UnaryFunction f) {
  static_assert(isSquare<Sparsity11>());
  static_assert(isSquare<Sparsity22>());
  static_assert(Sparsity21::num_cols == Sparsity11::num_cols);
  static_assert(Sparsity22::num_rows == Sparsity21::num_rows);
  constexpr auto num_rows =
      Sparsity11::num_rows + Sparsity22::num_rows + Sparsity33::num_rows;

  const auto tree = computeEliminationTreeBlocked3x3(
      sparsity11, sparsity21, sparsity22, sparsity31, sparsity32, sparsity33);
  std::array<std::size_t, num_rows> visitor;
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < Sparsity11::num_rows; ++i) {
    foreachAncestorInSubtree(sparsity11, tree, i, visitor, f);
  }
  for (std::size_t i = 0; i < Sparsity21::num_rows; ++i) {
    constexpr auto row_offset = Sparsity11::num_rows;
    constexpr auto col_offset = Sparsity21::num_cols;
    foreachAncestorInSubtree(sparsity21, tree, i, visitor, f,
                             row_offset);
    foreachAncestorInSubtree(sparsity22, tree, i, visitor, f,
                             row_offset, col_offset);
  }
  for (std::size_t i = 0; i < Sparsity31::num_rows; ++i) {
    constexpr auto row_offset = Sparsity11::num_rows + Sparsity22::num_rows;
    constexpr auto col_offset2 = Sparsity11::num_cols;
    constexpr auto col_offset3 = Sparsity11::num_cols + Sparsity22::num_cols;
    foreachAncestorInSubtree(sparsity31, tree, i, visitor, f,
                             row_offset);
    foreachAncestorInSubtree(sparsity32, tree, i, visitor, f,
                             row_offset, col_offset2);
    foreachAncestorInSubtree(sparsity33, tree, i, visitor, f,
                             row_offset, col_offset3);
  }
}

template <class Sparsity11, class Sparsity21, class Sparsity22,
          class UnaryFunction>
constexpr void foreachNonZeroWithFillBlocked(const Sparsity11& sparsity11,
                                             const Sparsity21& sparsity21,
                                             const Sparsity22& sparsity22,
                                             UnaryFunction f) {
  foreachNonZeroWithFillBlocked3x3(
      sparsity11, sparsity21, sparsity22,
      makeEmptySparsityCSR<0, Sparsity11::num_cols>(),
      makeEmptySparsityCSR<0, Sparsity22::num_cols>(),
      makeEmptySparsityCSR<0, 0>(), f);
}

}  // namespace ctldl
