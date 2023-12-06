#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

template <std::size_t nnz_, std::size_t num_rows_, std::size_t num_cols_>
struct Sparsity {
  static constexpr auto nnz = nnz_;
  static constexpr auto num_rows = num_rows_;
  static constexpr auto num_cols = num_cols_;

  std::array<Entry, nnz> entries;

  template <class Entries>
  constexpr explicit Sparsity(const Entries& entries_init) {
    fixInitIfZeroLengthArray(entries);
    std::transform(std::cbegin(entries_init), std::cend(entries_init),
                   std::begin(entries), [](const auto& entry) {
                     return Entry{entry.row_index, entry.col_index};
                   });
  }
};

template <std::size_t num_rows, std::size_t num_cols, class Entries>
constexpr auto makeSparsity(const Entries& entries) {
  using std::size;
  constexpr auto nnz = std::size_t{size(Entries{})};
  return Sparsity<nnz, num_rows, num_cols>(entries);
}

template <std::size_t num_rows, std::size_t num_cols, class Entry,
          std::size_t nnz>
constexpr auto makeSparsity(const Entry (&entries)[nnz]) {
  return Sparsity<nnz, num_rows, num_cols>(entries);
}

template <std::size_t num_rows, std::size_t num_cols, std::size_t nnz>
constexpr auto makeSparsity(const Entry (&entries)[nnz]) {
  return Sparsity<nnz, num_rows, num_cols>(entries);
}

template <class SparsityIn>
constexpr auto makeSparsity(const SparsityIn& sparsity) {
  return makeSparsity<SparsityIn::num_rows, SparsityIn::num_cols>(
      sparsity.entries);
}

template <std::size_t num_rows, std::size_t num_cols>
constexpr auto makeEmptySparsity() {
  return makeSparsity<num_rows, num_cols>(std::array<Entry, 0>{});
}

}  // namespace ctldl
