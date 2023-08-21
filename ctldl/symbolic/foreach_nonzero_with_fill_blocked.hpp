#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/compute_elimination_tree_blocked.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity11, class Sparsity21, class Sparsity22,
          class UnaryFunction>
constexpr void foreachNonZeroWithFillBlocked(const Sparsity11& sparsity11,
                                             const Sparsity21& sparsity21,
                                             const Sparsity22& sparsity22,
                                             UnaryFunction f) {
  static_assert(isSquare<Sparsity11>());
  static_assert(isSquare<Sparsity22>());
  static_assert(Sparsity21::num_cols == Sparsity11::num_cols);
  static_assert(Sparsity22::num_rows == Sparsity21::num_rows);
  constexpr auto num_rows = Sparsity11::num_rows + Sparsity21::num_rows;

  const auto tree =
      computeEliminationTreeBlocked(sparsity11, sparsity21, sparsity22);
  std::array<std::size_t, num_rows> visitor;
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < Sparsity11::num_rows; ++i) {
    foreachAncestorInSubtree(sparsity11, tree, i, visitor, f);
  }
  constexpr auto row_offset = Sparsity11::num_rows;
  constexpr auto col_offset = Sparsity21::num_cols;
  for (std::size_t i = 0; i < Sparsity21::num_rows; ++i) {
    foreachAncestorInSubtree(sparsity21, tree, i, visitor, f,
                             row_offset);
    foreachAncestorInSubtree(sparsity22, tree, i, visitor, f,
                             row_offset, col_offset);
  }
}

}  // namespace ctldl
