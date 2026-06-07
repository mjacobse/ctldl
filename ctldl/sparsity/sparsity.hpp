#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/utility/contracts.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>
#include <ctldl/utility/static_sized_range.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <meta>
#include <ranges>
#include <span>
#include <vector>

namespace ctldl {

template <std::size_t nnz_, std::size_t num_rows_, std::size_t num_cols_>
struct SparsityStatic {
  static constexpr std::size_t nnz() { return nnz_; }
  static constexpr std::size_t numRows() { return num_rows_; }
  static constexpr std::size_t numCols() { return num_cols_; }

  // only public to allow usage as NTTP
  std::array<Entry, nnz_> m_entries_do_not_touch;
  constexpr const std::array<Entry, nnz_>& entries() const {
    return m_entries_do_not_touch;
  }

  constexpr explicit SparsityStatic(
      std::ranges::input_range auto&& entries_init) {
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
    requires(
        std::ranges::input_range<decltype(SparsityIn::entries())> &&
        std::convertible_to<decltype(SparsityIn::numRows()), std::size_t> &&
        std::convertible_to<decltype(SparsityIn::numCols()), std::size_t>)
  constexpr SparsityStatic(const SparsityIn& sparsity_in)
      : SparsityStatic(sparsity_in.entries()) {}
};

template <class SparsityIn>
SparsityStatic(const SparsityIn& sparsity_in)
    -> SparsityStatic<std::tuple_size_v<decltype(SparsityIn::entries())>,
                      SparsityIn::numRows(), SparsityIn::numCols()>;

template <std::size_t num_rows, std::size_t num_cols,
          static_sized_range Entries>
constexpr auto makeSparsityStatic(const Entries& entries) {
  constexpr auto nnz = range_static_size_v<Entries>;
  return SparsityStatic<nnz, num_rows, num_cols>(entries);
}

template <std::size_t num_rows, std::size_t num_cols, class Entry,
          std::size_t nnz>
constexpr auto makeSparsityStatic(const Entry (&entries)[nnz]) {
  return SparsityStatic<nnz, num_rows, num_cols>(entries);
}

template <std::size_t num_rows, std::size_t num_cols, std::size_t nnz>
constexpr auto makeSparsityStatic(const Entry (&entries)[nnz]) {
  return SparsityStatic<nnz, num_rows, num_cols>(entries);
}

template <class SparsityIn>
constexpr auto makeSparsityStatic(const SparsityIn& sparsity) {
  return makeSparsityStatic<SparsityIn::numRows(), SparsityIn::numCols()>(
      sparsity.entries());
}

template <std::size_t num_rows, std::size_t num_cols>
constexpr auto makeEmptySparsityStatic() {
  return makeSparsityStatic<num_rows, num_cols>(std::array<Entry, 0>{});
}

class SparsityDynamic {
 private:
  static constexpr std::vector<Entry> makeEntries(
      std::ranges::input_range auto&& entries_init) {
    const auto entries_transformed =
        std::views::transform(entries_init, [](const auto& entry) {
          return Entry{entry.row_index, entry.col_index};
        });
    return std::vector<Entry>(entries_transformed.begin(),
                              entries_transformed.end());
  }

  std::vector<Entry> m_entries;
  std::size_t m_num_rows;
  std::size_t m_num_cols;

 public:
  constexpr const std::vector<Entry>& entries() const { return m_entries; }
  constexpr std::size_t numRows() const { return m_num_rows; }
  constexpr std::size_t numCols() const { return m_num_cols; }
  constexpr std::size_t nnz() const { return m_entries.size(); }

  constexpr explicit SparsityDynamic(
      const std::size_t num_rows_, const std::size_t num_cols_,
      std::ranges::input_range auto&& entries_init)
      : m_entries(makeEntries(entries_init)),
        m_num_rows(num_rows_),
        m_num_cols(num_cols_) {
    pre(std::ranges::all_of(
        entries_init, [num_rows_, num_cols_](const auto entry) {
          return entry.row_index < num_rows_ && entry.col_index < num_cols_;
        }));
  }

  template <std::size_t nnz_, std::size_t num_rows_, std::size_t num_cols_>
  explicit constexpr SparsityDynamic(
      const SparsityStatic<nnz_, num_rows_, num_cols_>& sparsity)
      : SparsityDynamic(sparsity.numRows(), sparsity.numCols(),
                        sparsity.entries()) {}
};

constexpr auto makeEmptySparsityDynamic(const std::size_t num_rows,
                                        const std::size_t num_cols) {
  return SparsityDynamic(num_rows, num_cols, std::vector<Entry>{});
}

class SparsityView {
 public:
  constexpr SparsityView(const SparsityDynamic& sparsity)
      : m_entries(sparsity.entries()),
        m_num_rows(sparsity.numRows()),
        m_num_cols(sparsity.numCols()) {}

  template <std::size_t nnz_, std::size_t num_rows_, std::size_t num_cols_>
  constexpr SparsityView(
      const SparsityStatic<nnz_, num_rows_, num_cols_>& sparsity)
      : m_entries(sparsity.entries()),
        m_num_rows(sparsity.numRows()),
        m_num_cols(sparsity.numCols()) {}

  constexpr const std::span<const Entry> entries() const { return m_entries; }
  constexpr std::size_t numRows() const { return m_num_rows; }
  constexpr std::size_t numCols() const { return m_num_cols; }
  constexpr std::size_t nnz() const { return m_entries.size(); }

 protected:
  constexpr SparsityView(const std::size_t num_rows, const std::size_t num_cols,
                         const std::span<const Entry> entries)
      : m_entries(entries), m_num_rows(num_rows), m_num_cols(num_cols) {
    pre(std::ranges::all_of(entries, [num_rows, num_cols](const auto entry) {
      return entry.row_index < num_rows && entry.col_index < num_cols;
    }));
  }

  friend consteval auto defineStaticSparsity(const SparsityView sparsity);

 private:
  std::span<const Entry> m_entries;
  std::size_t m_num_rows;
  std::size_t m_num_cols;
};

consteval auto defineStaticSparsity(const SparsityView sparsity) {
  const auto entries = std::define_static_array(sparsity.entries());
  return SparsityView(sparsity.numRows(), sparsity.numCols(), entries);
}

namespace detail {

template <std::size_t num_rows, std::size_t num_cols, std::size_t nnz,
          const Entry* data>
inline constexpr auto sparsity_static_helper_to_reflect_on =
    makeSparsityStatic<num_rows, num_cols>(
        std::span<const Entry, nnz>(data, nnz));

}  // namespace detail

consteval auto reflectSparsityStatic(const SparsityView sparsity) {
  const std::meta::info num_rows =
      std::meta::reflect_constant(sparsity.numRows());
  const std::meta::info num_cols =
      std::meta::reflect_constant(sparsity.numCols());
  const std::meta::info nnz = std::meta::reflect_constant(sparsity.nnz());
  const std::meta::info data =
      std::meta::reflect_constant(sparsity.entries().data());
  return std::meta::substitute(^^detail::sparsity_static_helper_to_reflect_on,
                               {num_rows, num_cols, nnz, data});
}

}  // namespace ctldl
