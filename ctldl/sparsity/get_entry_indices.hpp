#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_nonzero_info.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <std::size_t num_rows, std::size_t num_cols, std::size_t nnz>
constexpr auto getEntryIndices(const std::array<Entry, nnz> entries) {
  std::array<std::array<std::size_t, num_cols>, num_rows> entry_indices{};
  for (std::size_t entry_index = 0; entry_index < nnz; ++entry_index) {
    const std::size_t row_index = entries[entry_index].row_index;
    const std::size_t col_index = entries[entry_index].col_index;
    entry_indices[row_index][col_index] = entry_index;
  }
  return entry_indices;
}

}  // namespace ctldl
