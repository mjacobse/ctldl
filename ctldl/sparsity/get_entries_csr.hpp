#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/contracts.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

namespace ctldl {

constexpr auto getRowCounts(const SparsityView sparsity) {
  std::vector<std::size_t> row_counts(sparsity.numRows(), 0);
  for (const auto entry : sparsity.entries()) {
    row_counts[entry.row_index] += 1;
  }
  return row_counts;
}

constexpr auto getRowBeginIndices(
    const std::span<const std::size_t> row_counts) {
  pre(row_counts.size() < std::numeric_limits<std::size_t>::max());

  std::vector<std::size_t> row_begin_indices(row_counts.size() + 1);
  contract_assert(!row_begin_indices.empty());
  row_begin_indices[0] = 0;
  std::partial_sum(row_counts.begin(), row_counts.end(),
                   row_begin_indices.begin() + 1);
  return row_begin_indices;
}

struct EntriesCSR {
  std::vector<Entry> entries;
  std::vector<std::size_t> row_begin_indices;
};

constexpr auto getEntriesCSR(const SparsityView sparsity) {
  const auto row_counts = getRowCounts(sparsity);
  const auto row_begin_indices = getRowBeginIndices(row_counts);

  std::vector<Entry> entries(sparsity.nnz());
  auto row_insert_indices = row_begin_indices;
  for (const auto entry : sparsity.entries()) {
    entries[row_insert_indices[entry.row_index]] =
        Entry{entry.row_index, entry.col_index};
    row_insert_indices[entry.row_index] += 1;
  }
  return EntriesCSR{entries, row_begin_indices};
}

template <class Sparsity>
constexpr auto getEntriesStaticCSR(const Sparsity& sparsity) {
  std::array<Entry, Sparsity::nnz()> entries;
  std::array<std::size_t, Sparsity::numRows() + 1> row_begin_indices;
  fixInitIfZeroLengthArray(entries);
  fixInitIfZeroLengthArray(row_begin_indices);

  const auto entries_csr = getEntriesCSR(sparsity);
  std::ranges::copy(entries_csr.entries, entries.begin());
  std::ranges::copy(entries_csr.row_begin_indices, row_begin_indices.begin());
  return std::pair{entries, row_begin_indices};
}

}  // namespace ctldl
