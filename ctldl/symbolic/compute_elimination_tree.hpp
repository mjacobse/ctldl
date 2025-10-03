#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity>
constexpr auto computeEliminationTree(const Sparsity& sparsity) {
  static_assert(isSquare<Sparsity>());
  constexpr auto dim = Sparsity::numRows();

  EliminationTree<dim> tree;
  std::array<std::size_t, dim> ancestors;
  std::iota(ancestors.begin(), ancestors.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    addRowEliminationTree(sparsity, i, tree, ancestors);
  }
  return tree;
}

}  // namespace ctldl
