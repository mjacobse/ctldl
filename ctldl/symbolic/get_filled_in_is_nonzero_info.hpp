#pragma once

#include <ctldl/sparsity/is_nonzero_info.hpp>
#include <ctldl/symbolic/compute_elimination_tree.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity>
constexpr auto getFilledInIsNonzeroInfo() {
  static_assert(Sparsity::num_rows == Sparsity::num_cols);
  constexpr auto dim = Sparsity::num_rows;

  IsNonzeroInfo<dim, dim> is_nonzero;
  const auto mark_nonzero = [&is_nonzero](const std::size_t i,
                                          const std::size_t j) {
    is_nonzero[i][j] = true;
  };

  constexpr auto tree = computeEliminationTree<Sparsity>();
  std::array<std::size_t, dim> visitor;
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    foreachAncestorInSubtree<Sparsity>(tree, i, visitor, mark_nonzero);
  }

  return is_nonzero;
}

}  // namespace ctldl
