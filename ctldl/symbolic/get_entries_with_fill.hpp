#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>
#include <ctldl/symbolic/foreach_nonzero_with_fill.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class Sparsity>
constexpr auto getNumNonZerosWithFill(const Sparsity& sparsity) {
  std::size_t nnz = 0;
  const auto count_nonzero = [&nnz](const std::size_t /*i*/,
                                    const std::size_t /*j*/) { nnz += 1; };
  foreachNonZeroWithFill(sparsity, count_nonzero);
  return nnz;
}

template <auto sparsity>
constexpr auto getEntriesWithFill() {
  constexpr auto nnz = getNumNonZerosWithFill(sparsity);
  std::array<Entry, nnz> entries;
  fixInitIfZeroLengthArray(entries);
  std::size_t entry_index = 0;
  const auto add_entry = [&entries, &entry_index](const std::size_t i,
                                                  const std::size_t j) {
    entries[entry_index] = Entry{i, j};
    entry_index += 1;
  };
  foreachNonZeroWithFill(sparsity, add_entry);
  // sorting is not needed for correctness, but helps performance
  sortEntriesRowMajorSorted(entries);
  return entries;
}

}  // namespace ctldl
