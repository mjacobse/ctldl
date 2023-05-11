#pragma once

#include <ctldl/symbolic/compute_elimination_tree_repeating.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class SparsityA, class SparsityB, class UnaryFunction>
constexpr void foreachNonZeroWithFillRepeated(UnaryFunction f) {
  static_assert(SparsityA::num_rows == SparsityA::num_cols);
  static_assert(SparsityB::num_rows == SparsityA::num_rows);
  static_assert(SparsityB::num_cols == SparsityA::num_cols);
  constexpr auto dim = std::size_t{SparsityA::num_rows};

  constexpr auto tree = computeEliminationTreeRepeating<SparsityA, SparsityB>();
  std::array<std::size_t, 2 * dim> visitor;
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    constexpr auto row_offset = dim;
    constexpr auto col_offset = dim;
    foreachAncestorInSubtree<SparsityB>(tree, i, visitor, f, row_offset);
    foreachAncestorInSubtree<SparsityA>(tree, i, visitor, f, row_offset,
                                        col_offset);
  }
}

}  // namespace ctldl
