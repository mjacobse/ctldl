#pragma once

#include <ctldl/sparsity/get_entries_csr.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>

namespace ctldl {

template <std::size_t nnz, std::size_t num_rows, std::size_t num_cols>
struct SparsityCSR : public Sparsity<nnz, num_rows, num_cols> {
 private:
  using Base = Sparsity<nnz, num_rows, num_cols>;

  constexpr explicit SparsityCSR(
      std::pair<std::array<Entry, nnz>, std::array<std::size_t, num_rows + 1>>
          helper)
      : Base(helper.first), row_begin_indices(helper.second) {}

 public:
  std::array<std::size_t, num_rows + 1> row_begin_indices;

  template <class SparsityIn>
  constexpr explicit SparsityCSR(const SparsityIn& sparsity)
      : SparsityCSR(getEntriesCSR(sparsity)) {}

  constexpr auto find(const std::size_t i, const std::size_t j) const {
    const auto it_row_begin = Base::entries.cbegin() + row_begin_indices[i];
    const auto it_row_end = Base::entries.cbegin() + row_begin_indices[i + 1];
    return std::find_if(it_row_begin, it_row_end,
                        [j](const auto entry) { return entry.col_index == j; });
  }

  constexpr bool isNonZero(const std::size_t i, const std::size_t j) const {
    const auto it_row_end = Base::entries.cbegin() + row_begin_indices[i + 1];
    return find(i, j) != it_row_end;
  }

  constexpr auto entryIndex(const std::size_t i, const std::size_t j) const {
    assert(isNonZero(i, j));
    return static_cast<std::size_t>(find(i, j) - Base::entries.cbegin());
  }
};

template <class SparsityIn>
SparsityCSR(const SparsityIn&)
    -> SparsityCSR<SparsityIn::nnz, SparsityIn::num_rows, SparsityIn::num_cols>;

template <class SparsityIn>
constexpr auto makeSparsityCSR(const SparsityIn& sparsity) {
  return SparsityCSR{makeSparsity(sparsity)};
}

template <std::size_t num_rows, std::size_t num_cols, class Entries>
constexpr auto makeSparsityCSR(const Entries& entries) {
  return SparsityCSR{makeSparsity<num_rows, num_cols>(entries)};
}

template <std::size_t num_rows, std::size_t num_cols>
constexpr auto makeEmptySparsityCSR() {
  return SparsityCSR{makeEmptySparsity<num_rows, num_cols>()};
}

}  // namespace ctldl
