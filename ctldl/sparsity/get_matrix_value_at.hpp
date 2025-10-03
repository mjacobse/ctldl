#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>


namespace ctldl {

template <class Sparsity>
constexpr std::optional<std::size_t> findEntryIndex(const Sparsity& sparsity,
                                                    const Entry entry) {
  const auto& entries = sparsity.entries();
  const auto entries_begin = std::cbegin(entries);
  const auto entries_end = std::cend(entries);
  const auto it = std::find(entries_begin, entries_end, entry);
  if (it == entries_end) {
    return std::nullopt;
  }
  return static_cast<std::size_t>(it - entries_begin);
}

template <std::size_t i, std::size_t j, class Matrix>
constexpr auto getMatrixValueAt(const Matrix& matrix) {
  using Value = decltype(matrix.valueAt(0));

  constexpr std::optional<std::size_t> entry_index =
      findEntryIndex(Matrix::sparsity, Entry{i, j});
  if constexpr (entry_index.has_value()) {
    return matrix.valueAt(entry_index.value());
  }
  return Value{0.0};
}

}  // namespace ctldl
