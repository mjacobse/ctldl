#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/compute_elimination_tree.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>
#include <numeric>
#include <vector>

namespace ctldl {

template <class UnaryFunction>
constexpr void foreachNonZeroWithFill(const SparsityViewCSR sparsity,
                                      UnaryFunction f) {
  pre(isSquare(sparsity));
  const auto dim = sparsity.numRows();

  const auto tree = computeEliminationTree(sparsity);
  std::vector<std::size_t> visitor(dim);
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    foreachAncestorInSubtree(sparsity, tree, i, visitor, f);
  }
}

}  // namespace ctldl
