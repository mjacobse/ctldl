#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>

namespace ctldl {

template <class IsNonzeroGetter>
constexpr auto getEntries(IsNonzeroGetter getter) {
  constexpr auto is_nonzero = getter();
  constexpr auto num_rows = is_nonzero.num_rows;
  constexpr auto num_cols = is_nonzero.num_cols;
  constexpr int nnz = is_nonzero.nnz();
  std::array<Entry, nnz> entries{};
  int entry_index = 0;
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      if (is_nonzero[i][j]) {
        entries[entry_index] = Entry{i, j};
        entry_index += 1;
      }
    }
  }
  return entries;
}

template <class IsNonzeroGetter>
constexpr auto getEntriesCSR(IsNonzeroGetter getter) {
  constexpr auto is_nonzero = getter();
  constexpr auto num_rows = is_nonzero.num_rows;
  constexpr auto num_cols = is_nonzero.num_cols;
  constexpr int nnz = is_nonzero.nnz();
  std::array<Entry, nnz> entries{};
  std::array<int, num_rows + 1> row_begin_indices{};
  int entry_index = 0;
  for (int i = 0; i < num_rows; ++i) {
    row_begin_indices[i] = entry_index;
    for (int j = 0; j < num_cols; ++j) {
      if (is_nonzero[i][j]) {
        entries[entry_index] = Entry{i, j};
        entry_index += 1;
      }
    }
  }
  row_begin_indices.back() = entry_index;
  return std::pair{entries, row_begin_indices};
}

}  // namespace ctldl
