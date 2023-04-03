#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_nonzero_info.hpp>

#include <array>

namespace ctldl {

template <int num_rows, int num_cols, int nnz>
constexpr auto getEntryIndices(const std::array<Entry, nnz> entries) {
  std::array<std::array<int, num_cols>, num_rows> entry_indices{};
  for (int entry_index = 0; entry_index < nnz; ++entry_index) {
    const auto row_index = entries[entry_index].row_index;
    const auto col_index = entries[entry_index].col_index;
    entry_indices[row_index][col_index] = entry_index;
  }
  return entry_indices;
}

}  // namespace ctldl
