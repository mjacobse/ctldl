#pragma once

#include <ctldl/sparsity/get_entries_csr.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <span>

namespace ctldl {

template <std::size_t nnz, std::size_t num_rows, std::size_t num_cols>
struct SparsityCSR : public Sparsity<nnz, num_rows, num_cols> {
 private:
  using Base = Sparsity<nnz, num_rows, num_cols>;

  constexpr explicit SparsityCSR(
      std::pair<std::array<Entry, nnz>, std::array<std::size_t, num_rows + 1>>
          helper)
      : Base(helper.first), m_row_begin_indices(helper.second) {}

 public:
  std::array<std::size_t, num_rows + 1> m_row_begin_indices;
  constexpr const std::array<std::size_t, num_rows + 1>& rowBeginIndices()
      const {
    return m_row_begin_indices;
  }

  template <class SparsityIn>
  constexpr explicit SparsityCSR(const SparsityIn& sparsity)
      : SparsityCSR(getEntriesCSR(sparsity)) {}

  constexpr auto rowView(const std::size_t i) const {
    return std::span{Base::entries().data() + rowBeginIndices()[i],
                     Base::entries().data() + rowBeginIndices()[i + 1]};
  }

  constexpr auto find(const std::size_t i, const std::size_t j) const {
    const auto row = rowView(i);
    return std::find_if(std::cbegin(row), std::cend(row),
                        [j](const auto entry) { return entry.col_index == j; });
  }

  constexpr bool isNonZero(const std::size_t i, const std::size_t j) const {
    return find(i, j) != std::cend(rowView(i));
  }

  constexpr auto entryIndex(const std::size_t i, const std::size_t j) const {
    assert(isNonZero(i, j));
    return static_cast<std::size_t>(find(i, j) - std::cbegin(rowView(0)));
  }
};

template <class SparsityIn>
SparsityCSR(const SparsityIn&)
    -> SparsityCSR<SparsityIn::nnz(), SparsityIn::numRows(),
                   SparsityIn::numCols()>;

template <class SparsityIn>
constexpr auto makeSparsityCSR(const SparsityIn& sparsity) {
  return SparsityCSR{makeSparsity(sparsity)};
}

template <std::size_t num_rows, std::size_t num_cols, class Entries>
constexpr auto makeSparsityCSR(const Entries& entries) {
  return SparsityCSR{makeSparsity<num_rows, num_cols>(entries)};
}

template <std::size_t num_rows, std::size_t num_cols, class Entry,
          std::size_t nnz>
constexpr auto makeSparsityCSR(const Entry (&entries)[nnz]) {
  return SparsityCSR{makeSparsity<num_rows, num_cols>(entries)};
}

template <std::size_t num_rows, std::size_t num_cols, std::size_t nnz>
constexpr auto makeSparsityCSR(const Entry (&entries)[nnz]) {
  return SparsityCSR{makeSparsity<num_rows, num_cols>(entries)};
}

template <std::size_t num_rows, std::size_t num_cols>
constexpr auto makeEmptySparsityCSR() {
  return SparsityCSR{makeEmptySparsity<num_rows, num_cols>()};
}

}  // namespace ctldl
