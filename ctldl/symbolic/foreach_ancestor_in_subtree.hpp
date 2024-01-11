#pragma once

#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>


namespace ctldl {

template <class Sparsity, std::size_t dim, class GetParentFunction,
          class BinaryFunction>
constexpr void foreachAncestorInSubtreeImpl(
    const Sparsity& sparsity, const GetParentFunction get_parent,
    const std::size_t row_index, std::array<std::size_t, dim>& visitor,
    BinaryFunction f, const std::size_t row_offset = 0,
    const std::size_t col_offset = 0) {
  const auto i = row_index + row_offset;
  for (const auto entry : sparsity.rowView(row_index)) {
    auto j = entry.col_index + col_offset;
    while (visitor[j] != i) {
      f(i, j);
      visitor[j] = i;
      j = get_parent(j);
    }
  }
}

template <class Sparsity, std::size_t dim, class BinaryFunction>
constexpr void foreachAncestorInSubtree(const Sparsity& sparsity,
                                        const EliminationTree<dim>& tree,
                                        const std::size_t row_index,
                                        std::array<std::size_t, dim>& visitor,
                                        BinaryFunction f,
                                        const std::size_t row_offset = 0,
                                        const std::size_t col_offset = 0) {
  const auto get_parent = [&tree](const std::size_t j) {
    return tree.parent[j];
  };
  foreachAncestorInSubtreeImpl(sparsity, get_parent, row_index, visitor, f,
                               row_offset, col_offset);
}

template <class Sparsity, std::size_t dim, class BinaryFunction>
constexpr void foreachAncestorInSubtreeRepeating(
    const Sparsity& sparsity, const EliminationTree<2 * dim>& tree,
    const std::size_t row_index, std::array<std::size_t, dim>& visitor,
    BinaryFunction f, const std::size_t row_offset = 0,
    const std::size_t col_offset = 0) {
  const auto get_parent = [&tree](const std::size_t j) {
    return tree.parent[j] % dim;
  };
  foreachAncestorInSubtreeImpl(sparsity, get_parent, row_index, visitor, f,
                               row_offset, col_offset);
}

}  // namespace ctldl
