#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class IsNonzeroGetter>
constexpr auto getEntries(IsNonzeroGetter getter) {
  constexpr auto is_nonzero = getter();
  constexpr auto num_rows = std::size_t{is_nonzero.num_rows};
  constexpr auto num_cols = std::size_t{is_nonzero.num_cols};
  constexpr auto nnz = std::size_t{is_nonzero.nnz()};
  std::array<Entry, nnz> entries{};
  std::size_t entry_index = 0;
  for (std::size_t i = 0; i < num_rows; ++i) {
    for (std::size_t j = 0; j < num_cols; ++j) {
      if (is_nonzero[i][j]) {
        entries[entry_index] = Entry{i, j};
        entry_index += 1;
      }
    }
  }
  return entries;
}

template <class Sparsity>
constexpr auto getRowCounts() {
  std::array<std::size_t, Sparsity::num_rows> row_counts = {0};
  for (const auto entry : Sparsity::entries) {
    row_counts[entry.row_index] += 1;
  }
  return row_counts;
}

template <std::size_t dim>
constexpr auto getRowBeginIndices(
    const std::array<std::size_t, dim> row_counts) {
  std::array<std::size_t, dim + 1> row_begin_indices;
  row_begin_indices[0] = 0;
  std::partial_sum(row_counts.cbegin(), row_counts.cend(),
                   row_begin_indices.begin() + 1);
  return row_begin_indices;
}

template <class Sparsity>
constexpr auto getEntriesCSR() {
  constexpr auto row_counts = getRowCounts<Sparsity>();
  constexpr auto row_begin_indices = getRowBeginIndices(row_counts);
  constexpr auto nnz = std::size_t{std::size(Sparsity::entries)};

  std::array<Entry, nnz> entries;
  auto row_insert_indices = row_begin_indices;
  for (const auto entry : Sparsity::entries) {
    entries[row_insert_indices[entry.row_index]] = entry;
    row_insert_indices[entry.row_index] += 1;
  }
  return std::pair{entries, row_begin_indices};
}

}  // namespace ctldl
