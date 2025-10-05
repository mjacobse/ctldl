#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>
#include <numeric>
#include <vector>

namespace ctldl {

constexpr auto computeEliminationTree(const SparsityViewCSR sparsity) {
  pre(isSquare(sparsity));
  const auto dim = sparsity.numRows();

  EliminationTree tree(dim);
  std::vector<std::size_t> ancestors(dim);
  std::iota(ancestors.begin(), ancestors.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    addRowEliminationTree(sparsity, i, tree, ancestors);
  }
  return tree;
}

}  // namespace ctldl
