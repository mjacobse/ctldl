#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>
#include <numeric>


namespace ctldl {

template <class Sparsity11, class Sparsity21, class Sparsity22>
constexpr auto computeEliminationTreeBlocked(const Sparsity11& sparsity11,
                                             const Sparsity21& sparsity21,
                                             const Sparsity22& sparsity22) {
  static_assert(isSquare<Sparsity11>());
  static_assert(isSquare<Sparsity22>());
  static_assert(Sparsity21::num_cols == Sparsity11::num_cols);
  static_assert(Sparsity22::num_rows == Sparsity21::num_rows);
  constexpr auto num_rows = Sparsity11::num_rows + Sparsity21::num_rows;
  EliminationTree<num_rows> tree;

  std::array<std::size_t, num_rows> ancestors;
  std::iota(ancestors.begin(), ancestors.end(), 0);

  for (std::size_t i = 0; i < Sparsity11::num_rows; ++i) {
    addRowEliminationTree(sparsity11, i, tree, ancestors);
  }
  for (std::size_t i = 0; i < Sparsity21::num_rows; ++i) {
    constexpr auto row_offset = Sparsity11::num_rows;
    constexpr auto col_offset = Sparsity21::num_cols;
    addRowEliminationTree(sparsity21, i, tree, ancestors, row_offset);
    addRowEliminationTree(sparsity22, i, tree, ancestors, row_offset,
                          col_offset);
  }
  return tree;
}

}  // namespace ctldl
