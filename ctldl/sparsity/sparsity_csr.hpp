#pragma once

#include <ctldl/sparsity/get_entries_csr.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <span>

namespace ctldl {

template <std::size_t nnz, std::size_t num_rows, std::size_t num_cols>
struct SparsityStaticCSR : public SparsityStatic<nnz, num_rows, num_cols> {
 private:
  using Base = SparsityStatic<nnz, num_rows, num_cols>;

  constexpr explicit SparsityStaticCSR(
      std::pair<std::array<Entry, nnz>, std::array<std::size_t, num_rows + 1>>
          helper)
      : Base(helper.first), m_row_begin_indices_do_not_touch(helper.second) {}

 public:
  // only public to allow usage as NTTP
  std::array<std::size_t, num_rows + 1> m_row_begin_indices_do_not_touch;
  constexpr const std::array<std::size_t, num_rows + 1>& rowBeginIndices()
      const {
    return m_row_begin_indices_do_not_touch;
  }

  template <class SparsityIn>
  constexpr explicit SparsityStaticCSR(const SparsityIn& sparsity)
      : SparsityStaticCSR(getEntriesStaticCSR(sparsity)) {}

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
SparsityStaticCSR(const SparsityIn&)
    -> SparsityStaticCSR<SparsityIn::nnz(), SparsityIn::numRows(),
                         SparsityIn::numCols()>;

template <class SparsityIn>
constexpr auto makeSparsityStaticCSR(const SparsityIn& sparsity) {
  return SparsityStaticCSR{makeSparsityStatic(sparsity)};
}

template <std::size_t num_rows, std::size_t num_cols, class Entries>
constexpr auto makeSparsityStaticCSR(const Entries& entries) {
  return SparsityStaticCSR{makeSparsityStatic<num_rows, num_cols>(entries)};
}

template <std::size_t num_rows, std::size_t num_cols, class Entry,
          std::size_t nnz>
constexpr auto makeSparsityStaticCSR(const Entry (&entries)[nnz]) {
  return SparsityStaticCSR{makeSparsityStatic<num_rows, num_cols>(entries)};
}

template <std::size_t num_rows, std::size_t num_cols, std::size_t nnz>
constexpr auto makeSparsityStaticCSR(const Entry (&entries)[nnz]) {
  return SparsityStaticCSR{makeSparsityStatic<num_rows, num_cols>(entries)};
}

template <std::size_t num_rows, std::size_t num_cols>
constexpr auto makeEmptySparsityStaticCSR() {
  return SparsityStaticCSR{makeEmptySparsityStatic<num_rows, num_cols>()};
}

class SparsityDynamicCSR : public SparsityDynamic {
 private:
  std::vector<std::size_t> m_row_begin_indices = {0};

  constexpr explicit SparsityDynamicCSR(const std::size_t num_rows,
                                        const std::size_t num_cols,
                                        const EntriesCSR& entries_csr)
      : SparsityDynamic(num_rows, num_cols, entries_csr.entries),
        m_row_begin_indices(entries_csr.row_begin_indices) {}

 public:
  constexpr explicit SparsityDynamicCSR(const SparsityDynamic& sparsity)
      : SparsityDynamicCSR(sparsity.numRows(), sparsity.numCols(),
                           getEntriesCSR(sparsity)) {}

  template <std::size_t nnz, std::size_t num_rows_, std::size_t num_cols_>
  constexpr explicit SparsityDynamicCSR(
      const SparsityStaticCSR<nnz, num_rows_, num_cols_>& sparsity)
      : SparsityDynamic(sparsity.numRows(), sparsity.numCols(),
                        sparsity.entries()),
        m_row_begin_indices(sparsity.rowBeginIndices()) {}

  constexpr const std::vector<std::size_t>& rowBeginIndices() const {
    return m_row_begin_indices;
  }

  constexpr auto rowView(const std::size_t i) const {
    return std::span{entries().data() + m_row_begin_indices[i],
                     entries().data() + m_row_begin_indices[i + 1]};
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

constexpr auto makeEmptySparsityDynamicCSR(const std::size_t num_rows,
                                           const std::size_t num_cols) {
  return SparsityDynamicCSR(makeEmptySparsityDynamic(num_rows, num_cols));
}

class SparsityViewCSR : public SparsityView {
 private:
  std::span<const std::size_t> m_row_begin_indices;

 public:
  constexpr SparsityViewCSR(const SparsityDynamicCSR& sparsity)
      : SparsityView(sparsity),
        m_row_begin_indices(sparsity.rowBeginIndices()) {}

  template <std::size_t nnz, std::size_t num_rows_, std::size_t num_cols_>
  constexpr SparsityViewCSR(
      const SparsityStaticCSR<nnz, num_rows_, num_cols_>& sparsity)
      : SparsityView(sparsity),
        m_row_begin_indices(sparsity.rowBeginIndices()) {}

  constexpr const std::span<const std::size_t> rowBeginIndices() const {
    return m_row_begin_indices;
  }

  constexpr auto rowView(const std::size_t i) const {
    return std::span{entries().data() + m_row_begin_indices[i],
                     entries().data() + m_row_begin_indices[i + 1]};
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

}  // namespace ctldl
