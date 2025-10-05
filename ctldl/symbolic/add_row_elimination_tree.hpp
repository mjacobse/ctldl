#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>
#include <span>

namespace ctldl {

constexpr void addRowEliminationTree(const SparsityViewCSR sparsity,
                                     const std::size_t row_index,
                                     EliminationTree& tree,
                                     const std::span<std::size_t> ancestors,
                                     const std::size_t row_offset = 0,
                                     const std::size_t col_offset = 0) {
  pre(tree.parent.size() >= sparsity.numRows());

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
