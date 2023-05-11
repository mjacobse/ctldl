#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

template <std::size_t dim>
constexpr void sortEntriesRowMajorSorted(std::array<Entry, dim>& entries) {
  const auto row_major_sorted_order = [](const Entry lhs, const Entry rhs) {
    return lhs.row_index < rhs.row_index ||
           (lhs.row_index == rhs.row_index && lhs.col_index < rhs.col_index);
  };
  std::sort(entries.begin(), entries.end(), row_major_sorted_order);
}

}  // namespace ctldl
