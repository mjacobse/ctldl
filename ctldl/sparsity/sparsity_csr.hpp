#pragma once

#include <ctldl/sparsity/get_entries.hpp>
#include <ctldl/sparsity/get_entry_indices.hpp>
#include <ctldl/sparsity/get_is_nonzero_info.hpp>

#include <cassert>

namespace ctldl {

template <class Sparsity>
class SparsityCSR {
 public:
  static constexpr auto num_rows = Sparsity::num_rows;
  static constexpr auto num_cols = Sparsity::num_cols;
  static constexpr auto is_nonzero = getIsNonzeroInfo<Sparsity>();
 private:
  static constexpr auto m_helper = getEntriesCSR([] { return is_nonzero; });

 public:
  static constexpr auto entries = m_helper.first;
  static constexpr auto row_begin_indices = m_helper.second;
  static constexpr int nnz = entries.size();
 private:
  static constexpr auto m_entry_indices =
      getEntryIndices<num_rows, num_cols, nnz>(entries);

 public:
  static constexpr int entryIndex(const int i, const int j) {
    assert(is_nonzero[i][j]);
    assert(m_entry_indices[i][j] >= 0);
    return m_entry_indices[i][j];
  }
};

}  // namespace ctldl
