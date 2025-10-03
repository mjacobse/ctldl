#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/compute_elimination_tree_repeating.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity, std::size_t dim, class BinaryFunction>
constexpr void foreachAncestorInSubtreeRepeatingOuter(
    const Sparsity& sparsity, const EliminationTree<2 * dim>& tree,
    const std::size_t row_index, std::array<std::size_t, dim>& visitor,
    BinaryFunction f, const std::size_t row_offset = 0,
    const std::size_t col_offset = 0) {
  const auto get_parent = [&tree](const std::size_t j) {
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

template <class SparsityDiag, class SparsitySubdiag, class SparsityOuter,
          class UnaryFunction>
constexpr void foreachNonZeroWithFillRepeatingArrowhead(
    const SparsityDiag& sparsity_diag, const SparsitySubdiag& sparsity_subdiag,
    const SparsityOuter& sparsity_outer, UnaryFunction f) {
  static_assert(isSquare<SparsityDiag>());
  static_assert(isSquare<SparsitySubdiag>());
  static_assert(SparsitySubdiag::numRows() == SparsityDiag::numRows());
  static_assert(SparsityOuter::numCols() == SparsityDiag::numCols());
  constexpr auto dim = std::size_t{SparsityDiag::numRows()};
  constexpr auto dim_outer = std::size_t{SparsityOuter::numRows()};

  const auto tree =
      computeEliminationTreeRepeating(sparsity_diag, sparsity_subdiag);
  std::array<std::size_t, 2 * dim> visitor;
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < dim; ++i) {
    constexpr auto row_offset = dim;
    constexpr auto col_offset = dim;
    foreachAncestorInSubtree(sparsity_subdiag, tree, i, visitor, f, row_offset);
    foreachAncestorInSubtree(sparsity_diag, tree, i, visitor, f, row_offset,
                             col_offset);
  }

  std::array<std::size_t, dim> visitor_outer;
  std::iota(visitor_outer.begin(), visitor_outer.end(), 0);
  for (std::size_t i = 0; i < dim_outer; ++i) {
    constexpr auto row_offset = 2 * dim;
    foreachAncestorInSubtreeRepeatingOuter(sparsity_outer, tree, i,
                                           visitor_outer, f, row_offset);
  }
}

}  // namespace ctldl
