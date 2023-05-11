#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/symbolic/foreach_nonzero_with_fill.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity>
constexpr auto getNumNonZerosWithFill() {
  std::size_t nnz = 0;
  const auto count_nonzero = [&nnz](const std::size_t /*i*/,
                                    const std::size_t /*j*/) { nnz += 1; };
  foreachNonZeroWithFill<Sparsity>(count_nonzero);
  return nnz;
}

template <class Sparsity>
constexpr auto getEntriesWithFill() {
  constexpr auto nnz = getNumNonZerosWithFill<Sparsity>();
  std::array<Entry, nnz> entries;
  std::size_t entry_index = 0;
  const auto add_entry = [&entries, &entry_index](const std::size_t i,
                                                  const std::size_t j) {
    entries[entry_index] = Entry{i, j};
    entry_index += 1;
  };
  foreachNonZeroWithFill<Sparsity>(add_entry);
  return entries;
}

}  // namespace ctldl
