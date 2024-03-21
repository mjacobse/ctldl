#pragma once

#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity, std::size_t dim>
constexpr void addRowEliminationTree(const Sparsity& sparsity,
                                     const std::size_t row_index,
                                     EliminationTree<dim>& tree,
                                     std::array<std::size_t, dim>& ancestors,
                                     const std::size_t row_offset = 0,
                                     const std::size_t col_offset = 0) {
  static_assert(dim >= Sparsity::num_rows);

  const auto i = row_index + row_offset;
  for (const auto entry : sparsity.rowView(row_index)) {
    auto j = entry.col_index + col_offset;
    while (ancestors[j] < i) {
      if (!tree.hasParent(j)) {
        tree.parent[j] = i;
      }
      const auto current_ancestor = ancestors[j];
      ancestors[j] = i;
      j = current_ancestor;  // use shortcut to the top
    }
  }
}

}  // namespace ctldl
