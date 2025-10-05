#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/compute_elimination_tree_blocked.hpp>
#include <ctldl/symbolic/foreach_ancestor_in_subtree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class UnaryFunction>
constexpr void foreachNonZeroWithFillBlocked3x3(
    const SparsityViewCSR sparsity11, const SparsityViewCSR sparsity21,
    const SparsityViewCSR sparsity22, const SparsityViewCSR sparsity31,
    const SparsityViewCSR sparsity32, const SparsityViewCSR sparsity33,
    UnaryFunction f) {
  pre(isSquare(sparsity11));
  pre(isSquare(sparsity22));
  pre(sparsity21.numCols() == sparsity11.numCols());
  pre(sparsity22.numRows() == sparsity21.numRows());
  const auto num_rows =
      sparsity11.numRows() + sparsity22.numRows() + sparsity33.numRows();

  const auto tree = computeEliminationTreeBlocked3x3(
      sparsity11, sparsity21, sparsity22, sparsity31, sparsity32, sparsity33);
  std::vector<std::size_t> visitor(num_rows);
  std::iota(visitor.begin(), visitor.end(), 0);
  for (std::size_t i = 0; i < sparsity11.numRows(); ++i) {
    foreachAncestorInSubtree(sparsity11, tree, i, visitor, f);
  }
  for (std::size_t i = 0; i < sparsity21.numRows(); ++i) {
    const auto row_offset = sparsity11.numRows();
    const auto col_offset = sparsity21.numCols();
    foreachAncestorInSubtree(sparsity21, tree, i, visitor, f,
                             row_offset);
    foreachAncestorInSubtree(sparsity22, tree, i, visitor, f,
                             row_offset, col_offset);
  }
  for (std::size_t i = 0; i < sparsity31.numRows(); ++i) {
    const auto row_offset = sparsity11.numRows() + sparsity22.numRows();
    const auto col_offset2 = sparsity11.numCols();
    const auto col_offset3 = sparsity11.numCols() + sparsity22.numCols();
    foreachAncestorInSubtree(sparsity31, tree, i, visitor, f,
                             row_offset);
    foreachAncestorInSubtree(sparsity32, tree, i, visitor, f,
                             row_offset, col_offset2);
    foreachAncestorInSubtree(sparsity33, tree, i, visitor, f,
                             row_offset, col_offset3);
  }
}

template <class UnaryFunction>
constexpr void foreachNonZeroWithFillBlocked(const SparsityViewCSR sparsity11,
                                             const SparsityViewCSR sparsity21,
                                             const SparsityViewCSR sparsity22,
                                             UnaryFunction f) {
  foreachNonZeroWithFillBlocked3x3(
      sparsity11, sparsity21, sparsity22,
      makeEmptySparsityDynamicCSR(0, sparsity11.numCols()),
      makeEmptySparsityDynamicCSR(0, sparsity22.numCols()),
      makeEmptySparsityDynamicCSR(0, 0), f);
}

}  // namespace ctldl
