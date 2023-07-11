#pragma once

#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>


namespace ctldl {

template <class Sparsity, std::size_t dim, class BinaryFunction>
constexpr void foreachAncestorInSubtree(const Sparsity& sparsity,
                                        const EliminationTree<dim>& tree,
                                        const std::size_t i,
                                        std::array<std::size_t, dim>& visitor,
                                        BinaryFunction f,
                                        const std::size_t row_offset = 0,
                                        const std::size_t col_offset = 0) {
  for (const auto entry : sparsity.rowView(i)) {
    auto j = entry.col_index + col_offset;
    while (visitor[j] != i + row_offset) {
      f(i, j);
      visitor[j] = i + row_offset;
      j = tree.parent[j];
    }
  }
}

}  // namespace ctldl
