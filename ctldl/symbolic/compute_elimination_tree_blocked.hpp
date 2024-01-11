#pragma once

#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/add_row_elimination_tree.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>
#include <numeric>


namespace ctldl {

template <class Sparsity11, class Sparsity21, class Sparsity22,
          class Sparsity31, class Sparsity32, class Sparsity33>
constexpr auto computeEliminationTreeBlocked3x3(const Sparsity11& sparsity11,
                                                const Sparsity21& sparsity21,
                                                const Sparsity22& sparsity22,
                                                const Sparsity31& sparsity31,
                                                const Sparsity32& sparsity32,
                                                const Sparsity33& sparsity33) {
  static_assert(isSquare<Sparsity11>());
  static_assert(isSquare<Sparsity22>());
  static_assert(isSquare<Sparsity33>());
  static_assert(Sparsity21::num_cols == Sparsity11::num_cols);
  static_assert(Sparsity21::num_rows == Sparsity22::num_rows);
  static_assert(Sparsity31::num_cols == Sparsity11::num_cols);
  static_assert(Sparsity32::num_cols == Sparsity22::num_cols);
  static_assert(Sparsity31::num_rows == Sparsity33::num_rows);
  static_assert(Sparsity32::num_rows == Sparsity33::num_rows);
  constexpr auto num_rows =
      Sparsity11::num_rows + Sparsity22::num_rows + Sparsity33::num_rows;
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
  for (std::size_t i = 0; i < Sparsity31::num_rows; ++i) {
    constexpr auto row_offset = Sparsity11::num_rows + Sparsity22::num_rows;
    constexpr auto col_offset2 = Sparsity11::num_cols;
    constexpr auto col_offset3 = Sparsity11::num_cols + Sparsity22::num_cols;
    addRowEliminationTree(sparsity31, i, tree, ancestors, row_offset);
    addRowEliminationTree(sparsity32, i, tree, ancestors, row_offset,
                          col_offset2);
    addRowEliminationTree(sparsity33, i, tree, ancestors, row_offset,
                          col_offset3);
  }
  return tree;
}

template <class Sparsity11, class Sparsity21, class Sparsity22>
constexpr auto computeEliminationTreeBlocked(const Sparsity11& sparsity11,
                                             const Sparsity21& sparsity21,
                                             const Sparsity22& sparsity22) {
  return computeEliminationTreeBlocked3x3(
      sparsity11, sparsity21, sparsity22,
      makeEmptySparsityCSR<0, Sparsity11::num_cols>(),
      makeEmptySparsityCSR<0, Sparsity22::num_cols>(),
      makeEmptySparsityCSR<0, 0>());
}

}  // namespace ctldl
