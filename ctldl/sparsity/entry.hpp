#pragma once

#include <cstddef>

namespace ctldl {

struct Entry {
  std::size_t row_index;
  std::size_t col_index;
};

template <class OtherEntry>
constexpr bool operator==(const Entry& lhs, const OtherEntry& rhs) {
  // template parameter OtherEntry is so that we can also compare to any
  // user-type that has row_index and col_index
  return lhs.row_index == rhs.row_index && lhs.col_index == rhs.col_index;
}

}  // namespace ctldl
