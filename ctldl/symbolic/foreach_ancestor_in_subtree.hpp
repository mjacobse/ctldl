#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/elimination_tree.hpp>

#include <array>
#include <cstddef>
#include <span>


namespace ctldl {

template <class GetParentFunction, class BinaryFunction>
constexpr void foreachAncestorInSubtreeImpl(
    const SparsityViewCSR sparsity, const GetParentFunction get_parent,
    const std::size_t row_index, const std::span<std::size_t> visitor,
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

template <class BinaryFunction>
constexpr void foreachAncestorInSubtree(const SparsityViewCSR sparsity,
                                        const EliminationTree& tree,
                                        const std::size_t row_index,
                                        const std::span<std::size_t> visitor,
                                        BinaryFunction f,
                                        const std::size_t row_offset = 0,
                                        const std::size_t col_offset = 0) {
  const auto get_parent = [&tree](const std::size_t j) {
    return tree.getParent(j);
  };
  foreachAncestorInSubtreeImpl(sparsity, get_parent, row_index, visitor, f,
                               row_offset, col_offset);
}

}  // namespace ctldl
