#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

constexpr int getRepeatingMtxDim() { return 3; }

constexpr std::array<ctldl::Entry, 3> getRepeatingMtxEntriesA() {
  return {{{0, 0}, {1, 1}, {2, 2}}};
}

constexpr std::array<ctldl::Entry, 5> getRepeatingMtxEntriesB() {
  return {{{0, 0}, {1, 1}, {2, 1}, {1, 2}, {2, 2}}};
}

constexpr std::array<std::size_t, getRepeatingMtxDim()>
getRepeatingMtxPermutation() {
  return {0, 1, 2};
}
