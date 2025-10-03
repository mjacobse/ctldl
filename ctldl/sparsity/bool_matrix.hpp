#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>

namespace ctldl {

/**
 * A dense matrix containing boolean entries.
 *
 * Can be used to represent sparse matrices in a very inefficient, but sometimes
 * convenient way, where an entry is true if and only if there is a nonzero
 * entry in the matrix to be represented.
 */
template <std::size_t num_rows, std::size_t num_cols>
struct BoolMatrix {
  static constexpr std::size_t numRows() { return num_rows; }
  static constexpr std::size_t numCols() { return num_cols; }
  std::array<std::array<bool, numCols()>, numRows()> values;
  constexpr auto nnz() const {
    return std::accumulate(
        values.cbegin(), values.cend(), std::size_t{0},
        [](const std::size_t carry, const std::array<bool, numCols()>& row) {
          return carry +
                 static_cast<std::size_t>(std::ranges::count(row, true));
        });
  }
};

template <auto bool_matrix>
constexpr auto makeSparseEntriesFromBoolMatrix() {
  constexpr auto num_rows = std::size_t{bool_matrix.numRows()};
  constexpr auto num_cols = std::size_t{bool_matrix.numCols()};
  constexpr auto nnz = bool_matrix.nnz();

  std::array<Entry, nnz> entries;
  fixInitIfZeroLengthArray(entries);
  auto entry_index = std::size_t{0};
  for (std::size_t i = 0; i < num_rows; ++i) {
    for (std::size_t j = 0; j < num_cols; ++j) {
      if (bool_matrix.values[i][j]) {
        entries[entry_index] = Entry{i, j};
        entry_index += 1;
      }
    }
  }
  assert(entry_index == nnz);
  return entries;
}

}  // namespace ctldl
