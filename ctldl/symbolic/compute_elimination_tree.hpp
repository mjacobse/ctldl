#pragma once

#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity>
constexpr auto computeEliminationTree(const Sparsity& sparsity) {
  static_assert(Sparsity::num_rows == Sparsity::num_cols);
  constexpr auto dim = Sparsity::num_rows;

  EliminationTree<dim> tree;
  std::array<std::size_t, dim> ancestors;
  std::iota(ancestors.begin(), ancestors.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    addRowElimintationTree(sparsity, i, tree, ancestors);
  }
  return tree;
}

}  // namespace ctldl
