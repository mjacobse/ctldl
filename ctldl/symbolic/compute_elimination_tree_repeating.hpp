#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree_subtree.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <ranges>


namespace ctldl {

template <class SparsityA, class SparsityB>
constexpr auto computeEliminationTreeRepeating(const SparsityA& sparsity_A,
                                               const SparsityB& sparsity_B) {
  static_assert(isSquare<SparsityA>());
  static_assert(isSquare<SparsityB>());
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
      addRowEliminationTree(sparsity_B, i, tree, ancestors, row_offset);
      addRowEliminationTree(sparsity_A, i, tree, ancestors, row_offset,
                            col_offset);
    }

    tree_init = subtreeBack<dim>(tree);
    if (tree_init == tree_init_previous) {
      break;
    }

    setSubtreeFront(tree, tree_init);
    clearSubtreeBack(tree, dim);
    std::ranges::transform(
        std::views::drop(ancestors, dim), ancestors.begin(),
        [](const std::size_t ancestor) { return ancestor - dim; });
    std::iota(ancestors.begin() + dim, ancestors.end(), dim);
  };

  return tree;
}

}  // namespace ctldl
