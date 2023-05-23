#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class InputEntry, std::size_t num_rows, std::size_t num_cols>
constexpr Entry permutedEntry(const InputEntry entry,
                              const Permutation<num_rows>& permutation_row,
                              const Permutation<num_cols>& permutation_col) {
  return Entry{permutation_row[entry.row_index],
               permutation_col[entry.col_index]};
}

}  // namespace ctldl
