#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

template <class InputEntry, std::size_t dim>
constexpr Entry permutedEntryLowerTriangle(
    const InputEntry entry, const Permutation<dim>& permutation) {
  return Entry{
      std::max(permutation[entry.row_index], permutation[entry.col_index]),
      std::min(permutation[entry.row_index], permutation[entry.col_index])};
}

}  // namespace ctldl
