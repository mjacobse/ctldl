#pragma once

#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity, std::size_t dim>
constexpr void addRowElimintationTree(const std::size_t row_index,
                                      EliminationTree<dim>& tree,
                                      std::array<std::size_t, dim>& ancestors,
                                      const std::size_t row_offset = 0,
                                      const std::size_t col_offset = 0) {
  static_assert(dim >= Sparsity::num_rows);

  const auto row_begin = Sparsity::row_begin_indices[row_index];
  const auto row_end = Sparsity::row_begin_indices[row_index + 1];
  const auto i = row_index + row_offset;

  for (auto entry_index = row_begin; entry_index != row_end; ++entry_index) {
    auto j = Sparsity::entries[entry_index].col_index + col_offset;
    while (ancestors[j] < i) {
      if (!tree.hasParent(j)) {
        tree.parent[j] = i;
      }
      ancestors[j] = i;
      j = tree.parent[j];
    }
  }
}

}  // namespace ctldl
