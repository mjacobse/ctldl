#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <algorithm>
#include <span>

namespace ctldl {

constexpr void sortEntriesRowMajorSorted(const std::span<Entry> entries) {
  const auto row_major_sorted_order = [](const Entry lhs, const Entry rhs) {
    return lhs.row_index < rhs.row_index ||
           (lhs.row_index == rhs.row_index && lhs.col_index < rhs.col_index);
  };
  std::ranges::sort(entries, row_major_sorted_order);
}

}  // namespace ctldl
