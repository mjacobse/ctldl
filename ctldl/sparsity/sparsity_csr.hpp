#pragma once

#include <ctldl/sparsity/get_entries_csr.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>

namespace ctldl {

template <class Sparsity>
class SparsityCSR {
 private:
  static constexpr auto m_helper = getEntriesCSR<Sparsity>();
 public:
  static constexpr auto entries = m_helper.first;
  static constexpr auto row_begin_indices = m_helper.second;

  static constexpr auto num_rows = std::size_t{Sparsity::num_rows};
  static constexpr auto num_cols = std::size_t{Sparsity::num_cols};
  static constexpr auto nnz = std::size_t{entries.size()};

  static constexpr auto find(const std::size_t i, const std::size_t j) {
    const auto it_row_begin = entries.cbegin() + row_begin_indices[i];
    const auto it_row_end = entries.cbegin() + row_begin_indices[i + 1];
    return std::find_if(it_row_begin, it_row_end,
                        [j](const auto entry) { return entry.col_index == j; });
  }

  static constexpr bool isNonZero(const std::size_t i, const std::size_t j) {
    const auto it_row_end = entries.cbegin() + row_begin_indices[i + 1];
    return find(i, j) != it_row_end;
  }

  static constexpr auto entryIndex(const std::size_t i, const std::size_t j) {
    assert(isNonZero(i, j));
    return static_cast<std::size_t>(find(i, j) - entries.cbegin());
  }
};

}  // namespace ctldl
