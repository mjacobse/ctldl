#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <ranges>

namespace ctldl {

template <std::size_t nnz_, std::size_t num_rows_, std::size_t num_cols_>
struct Sparsity {
  static constexpr std::size_t nnz() { return nnz_; }
  static constexpr std::size_t numRows() { return num_rows_; }
  static constexpr std::size_t numCols() { return num_cols_; }

  // only public to allow usage as NTTP
  std::array<Entry, nnz_> m_entries_do_not_touch;
  constexpr const std::array<Entry, nnz_>& entries() const {
    return m_entries_do_not_touch;
  }

  constexpr explicit Sparsity(std::ranges::input_range auto&& entries_init) {
    fixInitIfZeroLengthArray(m_entries_do_not_touch);
    std::ranges::transform(entries_init, std::begin(m_entries_do_not_touch),
                           [](const auto& entry) {
                             return Entry{entry.row_index, entry.col_index};
                           });
    assert(std::ranges::all_of(entries(), [](const auto entry) {
      return entry.row_index < numRows() && entry.col_index < numCols();
    }));
  }

  template <class SparsityIn>
    requires(std::ranges::input_range<decltype(SparsityIn::entries())> &&
             std::convertible_to<decltype(SparsityIn::numRows()), std::size_t> &&
             std::convertible_to<decltype(SparsityIn::numCols()), std::size_t>)
  constexpr Sparsity(const SparsityIn& sparsity_in)
      : Sparsity(sparsity_in.entries()) {}
};

template <class SparsityIn>
Sparsity(const SparsityIn& sparsity_in)
    -> Sparsity<std::tuple_size_v<decltype(SparsityIn::entries())>,
                SparsityIn::numRows(), SparsityIn::numCols()>;

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
  return makeSparsity<SparsityIn::numRows(), SparsityIn::numCols()>(
      sparsity.entries());
}

template <std::size_t num_rows, std::size_t num_cols>
constexpr auto makeEmptySparsity() {
  return makeSparsity<num_rows, num_cols>(std::array<Entry, 0>{});
}

}  // namespace ctldl
