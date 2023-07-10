#pragma once

#include <ctldl/symbolic/compute_elimination_tree.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity, class UnaryFunction>
constexpr void foreachNonZeroWithFill(const Sparsity& sparsity,
                                      UnaryFunction f) {
  static_assert(Sparsity::num_rows == Sparsity::num_cols);
  constexpr auto dim = Sparsity::num_rows;

  const auto tree = computeEliminationTree(sparsity);
  std::array<std::size_t, dim> visitor;
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    foreachAncestorInSubtree(sparsity, tree, i, visitor, f);
  }
}

}  // namespace ctldl
