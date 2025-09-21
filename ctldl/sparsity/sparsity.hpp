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
  static constexpr auto nnz = nnz_;
  static constexpr auto num_rows = num_rows_;
  static constexpr auto num_cols = num_cols_;

  std::array<Entry, nnz> entries;

  constexpr explicit Sparsity(std::ranges::input_range auto&& entries_init) {
    fixInitIfZeroLengthArray(entries);
    std::transform(std::cbegin(entries_init), std::cend(entries_init),
                   std::begin(entries), [](const auto& entry) {
                     return Entry{entry.row_index, entry.col_index};
                   });
    assert(std::ranges::all_of(entries, [](const auto entry) {
      return entry.row_index < num_rows && entry.col_index < num_cols;
    }));
  }

  template <class SparsityIn>
    requires(std::ranges::input_range<decltype(SparsityIn::entries)> &&
             std::convertible_to<decltype(SparsityIn::num_rows), std::size_t> &&
             std::convertible_to<decltype(SparsityIn::num_cols), std::size_t>)
  constexpr Sparsity(const SparsityIn& sparsity_in)
      : Sparsity(sparsity_in.entries) {}
};

template <class SparsityIn>
Sparsity(const SparsityIn& sparsity_in)
    -> Sparsity<std::tuple_size_v<decltype(SparsityIn::entries)>,
                SparsityIn::num_rows, SparsityIn::num_cols>;

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
