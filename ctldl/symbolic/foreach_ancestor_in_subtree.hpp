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
  const auto row_begin = sparsity.row_begin_indices[i];
  const auto row_end = sparsity.row_begin_indices[i + 1];
  for (auto entry_index = row_begin; entry_index != row_end; ++entry_index) {
    auto j = sparsity.entries[entry_index].col_index + col_offset;
    while (visitor[j] != i + row_offset) {
      f(i, j);
      visitor[j] = i + row_offset;
      j = tree.parent[j];
    }
  }
}

}  // namespace ctldl
