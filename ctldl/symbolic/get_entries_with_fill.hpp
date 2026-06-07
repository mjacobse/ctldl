#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>
#include <ctldl/symbolic/foreach_nonzero_with_fill.hpp>

#include <cstddef>
#include <vector>

namespace ctldl {

constexpr auto getEntriesWithFill(const SparsityViewCSR sparsity) {
  std::vector<Entry> entries;
  const auto add_entry = [&entries](const std::size_t i, const std::size_t j) {
    entries.push_back(Entry{i, j});
  };
  foreachNonZeroWithFill(sparsity, add_entry);
  // sorting is not needed for correctness, but helps performance
  sortEntriesRowMajorSorted(entries);
  return entries;
}

}  // namespace ctldl
