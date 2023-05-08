#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>


namespace ctldl {

template <std::size_t i, std::size_t j, class Matrix>
constexpr auto getMatrixValueAt(const Matrix& matrix) {
  using Sparsity = typename Matrix::Sparsity;
  using Value = decltype(matrix.valueAt(0));

  constexpr auto entries_begin = std::cbegin(Sparsity::entries);
  constexpr auto entries_end = std::cend(Sparsity::entries);

  constexpr auto it =
      std::find_if(entries_begin, entries_end, [](const auto entry) {
        return entry.row_index == i && entry.col_index == j;
      });

  if constexpr (it != entries_end) {
    constexpr auto entry_index = std::size_t{it - entries_begin};
    return matrix.valueAt(entry_index);
  }
  return Value{0.0};
}

}  // namespace ctldl
