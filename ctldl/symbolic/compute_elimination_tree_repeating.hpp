#pragma once

#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree_subtree.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>


namespace ctldl {

template <class SparsityA, class SparsityB>
constexpr auto computeEliminationTreeRepeating() {
  static_assert(SparsityA::num_rows == SparsityA::num_cols);
  static_assert(SparsityB::num_rows == SparsityB::num_cols);
  static_assert(SparsityA::num_rows == SparsityB::num_rows);
  constexpr auto dim = SparsityA::num_rows;

  EliminationTree<2 * dim> tree;
  EliminationTree<dim> tree_init;

  std::array<std::size_t, 2 * dim> ancestors;
  std::iota(ancestors.begin(), ancestors.end(), 0);

  while (true) {
    const auto tree_init_previous = tree_init;
    for (std::size_t i = 0; i < dim; ++i) {
      constexpr auto row_offset = dim;
      constexpr auto col_offset = dim;
      addRowElimintationTree<SparsityB>(i, tree, ancestors, row_offset);
      addRowElimintationTree<SparsityA>(i, tree, ancestors, row_offset, col_offset);
    }

    tree_init = subtreeBack<dim>(tree);
    if (tree_init == tree_init_previous) {
      break;
    }

    setSubtreeFront(tree, tree_init);
    clearSubtreeBack(tree, dim);
    std::transform(ancestors.cbegin() + dim, ancestors.cend(),
                   ancestors.begin(),
                   [](const std::size_t ancestor) { return ancestor - dim; });
    std::iota(ancestors.begin() + dim, ancestors.end(), dim);
  };

  return tree;
}

}  // namespace ctldl
