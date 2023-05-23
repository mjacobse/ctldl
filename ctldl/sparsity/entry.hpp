#pragma once

#include <cstddef>

namespace ctldl {

struct Entry {
  std::size_t row_index;
  std::size_t col_index;
};

template <class EntryLhs, class EntryRhs>
constexpr bool operator==(const EntryLhs& lhs, const EntryRhs& rhs) {
  return lhs.row_index == rhs.row_index && lhs.col_index == rhs.col_index;
}

}  // namespace ctldl
