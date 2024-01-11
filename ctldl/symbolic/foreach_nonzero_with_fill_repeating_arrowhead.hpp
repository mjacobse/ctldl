#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/compute_elimination_tree_repeating.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class SparsityDiag, class SparsitySubdiag, class SparsityOuter,
          class UnaryFunction>
constexpr void foreachNonZeroWithFillRepeatingArrowhead(
    const SparsityDiag& sparsity_diag, const SparsitySubdiag& sparsity_subdiag,
    const SparsityOuter& sparsity_outer, UnaryFunction f) {
  static_assert(isSquare<SparsityDiag>());
  static_assert(isSquare<SparsitySubdiag>());
  static_assert(SparsitySubdiag::num_rows == SparsityDiag::num_rows);
  static_assert(SparsityOuter::num_cols == SparsityDiag::num_cols);
  constexpr auto dim = std::size_t{SparsityDiag::num_rows};
  constexpr auto dim_outer = std::size_t{SparsityOuter::num_rows};

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
    foreachAncestorInSubtreeRepeating(sparsity_outer, tree, i, visitor_outer, f,
                                      row_offset);
  }
}

}  // namespace ctldl
