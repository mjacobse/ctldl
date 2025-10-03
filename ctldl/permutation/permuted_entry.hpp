#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class InputEntry>
constexpr Entry permutedEntry(const InputEntry entry,
                              const PermutationView permutation_row,
                              const PermutationView permutation_col) {
  return Entry{permutation_row[entry.row_index],
               permutation_col[entry.col_index]};
}

}  // namespace ctldl
