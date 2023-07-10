#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <class Sparsity>
constexpr auto getRowCounts(const Sparsity& sparsity) {
  std::array<std::size_t, Sparsity::num_rows> row_counts;
  row_counts.fill(0);
  for (const auto entry : sparsity.entries) {
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
constexpr auto getEntriesCSR(const Sparsity& sparsity) {
  const auto row_counts = getRowCounts(sparsity);
  const auto row_begin_indices = getRowBeginIndices(row_counts);

  std::array<Entry, Sparsity::nnz> entries;
  auto row_insert_indices = row_begin_indices;
  for (const auto entry : sparsity.entries) {
    entries[row_insert_indices[entry.row_index]] =
        Entry{entry.row_index, entry.col_index};
    row_insert_indices[entry.row_index] += 1;
  }
  return std::pair{entries, row_begin_indices};
}

}  // namespace ctldl
