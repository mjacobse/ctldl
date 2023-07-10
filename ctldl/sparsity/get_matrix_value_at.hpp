#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>


namespace ctldl {

template <std::size_t i, std::size_t j, class Matrix>
constexpr auto getMatrixValueAt(const Matrix& matrix) {
  using Value = decltype(matrix.valueAt(0));

  constexpr auto entries_begin = std::cbegin(Matrix::sparsity.entries);
  constexpr auto entries_end = std::cend(Matrix::sparsity.entries);

  constexpr auto it = std::find(entries_begin, entries_end, Entry{i, j});
  if constexpr (it != entries_end) {
    constexpr auto entry_index = std::size_t{it - entries_begin};
    return matrix.valueAt(entry_index);
  }
  return Value{0.0};
}

}  // namespace ctldl
