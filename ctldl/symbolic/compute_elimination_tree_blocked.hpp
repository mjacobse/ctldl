#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>
#include <numeric>
#include <vector>


namespace ctldl {

constexpr auto computeEliminationTreeBlocked3x3(
    const SparsityViewCSR sparsity11, const SparsityViewCSR sparsity21,
    const SparsityViewCSR sparsity22, const SparsityViewCSR sparsity31,
    const SparsityViewCSR sparsity32, const SparsityViewCSR sparsity33) {
  pre(isSquare(sparsity11));
  pre(isSquare(sparsity22));
  pre(isSquare(sparsity33));
  pre(sparsity21.numCols() == sparsity11.numCols());
  pre(sparsity21.numRows() == sparsity22.numRows());
  pre(sparsity31.numCols() == sparsity11.numCols());
  pre(sparsity32.numCols() == sparsity22.numCols());
  pre(sparsity31.numRows() == sparsity33.numRows());
  pre(sparsity32.numRows() == sparsity33.numRows());
  const auto num_rows =
      sparsity11.numRows() + sparsity22.numRows() + sparsity33.numRows();
  EliminationTree tree(num_rows);

  std::vector<std::size_t> ancestors(num_rows);
  std::iota(ancestors.begin(), ancestors.end(), 0);

  for (std::size_t i = 0; i < sparsity11.numRows(); ++i) {
    addRowEliminationTree(sparsity11, i, tree, ancestors);
  }
  for (std::size_t i = 0; i < sparsity21.numRows(); ++i) {
    const auto row_offset = sparsity11.numRows();
    const auto col_offset = sparsity21.numCols();
    addRowEliminationTree(sparsity21, i, tree, ancestors, row_offset);
    addRowEliminationTree(sparsity22, i, tree, ancestors, row_offset,
                          col_offset);
  }
  for (std::size_t i = 0; i < sparsity31.numRows(); ++i) {
    const auto row_offset = sparsity11.numRows() + sparsity22.numRows();
    const auto col_offset2 = sparsity11.numCols();
    const auto col_offset3 = sparsity11.numCols() + sparsity22.numCols();
    addRowEliminationTree(sparsity31, i, tree, ancestors, row_offset);
    addRowEliminationTree(sparsity32, i, tree, ancestors, row_offset,
                          col_offset2);
    addRowEliminationTree(sparsity33, i, tree, ancestors, row_offset,
                          col_offset3);
  }
  return tree;
}

constexpr auto computeEliminationTreeBlocked(const SparsityViewCSR sparsity11,
                                             const SparsityViewCSR sparsity21,
                                             const SparsityViewCSR sparsity22) {
  return computeEliminationTreeBlocked3x3(
      sparsity11, sparsity21, sparsity22,
      makeEmptySparsityDynamicCSR(0, sparsity11.numCols()),
      makeEmptySparsityDynamicCSR(0, sparsity22.numCols()),
      makeEmptySparsityDynamicCSR(0, 0));
}

}  // namespace ctldl
