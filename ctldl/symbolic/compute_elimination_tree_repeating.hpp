#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree_subtree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <ranges>
#include <vector>


namespace ctldl {

constexpr auto computeEliminationTreeRepeating(
    const SparsityViewCSR sparsity_A, const SparsityViewCSR sparsity_B) {
  pre(isSquare(sparsity_A));
  pre(isSquare(sparsity_B));
  pre(sparsity_A.numRows() == sparsity_B.numRows());
  const auto dim = sparsity_A.numRows();

  EliminationTree tree(2 * dim);
  EliminationTree tree_init(dim);

  std::vector<std::size_t> ancestors(2 * dim);
  std::iota(ancestors.begin(), ancestors.end(), 0);

  while (true) {
    const auto tree_init_previous = tree_init;
    for (std::size_t i = 0; i < dim; ++i) {
      const auto row_offset = dim;
      const auto col_offset = dim;
      addRowEliminationTree(sparsity_B, i, tree, ancestors, row_offset);
      addRowEliminationTree(sparsity_A, i, tree, ancestors, row_offset,
                            col_offset);
    }

    tree_init = subtreeBack(tree, dim);
    if (tree_init == tree_init_previous) {
      break;
    }

    setSubtreeFront(tree, tree_init);
    clearSubtreeBack(tree, dim);
    std::ranges::transform(
        std::views::drop(ancestors, static_cast<std::ptrdiff_t>(dim)),
        ancestors.begin(),
        [dim](const std::size_t ancestor) { return ancestor - dim; });
    std::iota(ancestors.begin() + static_cast<std::ptrdiff_t>(dim),
              ancestors.end(), dim);
  };

  return tree;
}

}  // namespace ctldl
