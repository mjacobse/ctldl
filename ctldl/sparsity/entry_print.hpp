#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <ostream>

namespace ctldl {

inline std::ostream& operator<<(std::ostream& os, const Entry& right) {
  os << "Entry{.row_index=" << right.row_index
     << ", .col_index=" << right.col_index << '}';
  return os;
}

}  // namespace ctldl
