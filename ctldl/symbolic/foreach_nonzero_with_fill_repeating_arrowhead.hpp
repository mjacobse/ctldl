#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/compute_elimination_tree_repeating.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <array>
#include <cstddef>
#include <numeric>
#include <span>

namespace ctldl {

template <class BinaryFunction>
constexpr void foreachAncestorInSubtreeRepeatingOuter(
    const SparsityViewCSR sparsity, const EliminationTree& tree,
    const std::size_t row_index, const std::span<std::size_t> visitor,
    BinaryFunction f, const std::size_t row_offset = 0,
    const std::size_t col_offset = 0) {
  const auto dim = std::size_t{visitor.size()};
  const auto get_parent = [dim, &tree](const std::size_t j) {
    if (!tree.hasParent(j)) {
      // We use the elimination tree of the tridiagonal part (i.e. it does not
      // know anything about the outer part). So if the parent of the full tree
      // would be an outer node, then we will encounter no parent at all
      // instead. Since that would only lead to entries in the non-repeating
      // diagonal outer part which we do not care about right now, we can safely
      // skip them. To do that we can just return the node itself, which will
      // directly stop the iteration.
      return j;
    }
    // otherwise we want to treat columns of the subdiagonal block the same as
    // columns of the diagonal block, so we modulo them to the range of
    // [0, blocksize_tridiag).
    return tree.getParent(j) % dim;
  };
  foreachAncestorInSubtreeImpl(sparsity, get_parent, row_index, visitor, f,
                               row_offset, col_offset);
}

template <class UnaryFunction>
constexpr void foreachNonZeroWithFillRepeatingArrowhead(
    const SparsityViewCSR sparsity_diag, const SparsityViewCSR sparsity_subdiag,
    const SparsityViewCSR sparsity_outer, UnaryFunction f) {
  pre(isSquare(sparsity_diag));
  pre(isSquare(sparsity_subdiag));
  pre(sparsity_subdiag.numRows() == sparsity_diag.numRows());
  pre(sparsity_outer.numCols() == sparsity_diag.numCols());
  const auto dim = std::size_t{sparsity_diag.numRows()};
  const auto dim_outer = std::size_t{sparsity_outer.numRows()};

  const auto tree =
      computeEliminationTreeRepeating(sparsity_diag, sparsity_subdiag);
  std::vector<std::size_t> visitor(2 * dim);
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    const auto row_offset = dim;
    const auto col_offset = dim;
    foreachAncestorInSubtree(sparsity_subdiag, tree, i, visitor, f, row_offset);
    foreachAncestorInSubtree(sparsity_diag, tree, i, visitor, f, row_offset,
                             col_offset);
  }

  std::vector<std::size_t> visitor_outer(dim);
  std::iota(visitor_outer.begin(), visitor_outer.end(), 0);
  for (std::size_t i = 0; i < dim_outer; ++i) {
    const auto row_offset = 2 * dim;
    foreachAncestorInSubtreeRepeatingOuter(sparsity_outer, tree, i,
                                           visitor_outer, f, row_offset);
  }
}

}  // namespace ctldl
