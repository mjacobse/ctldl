#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

constexpr int getRepeatingMtxDim() { return 7; }

constexpr std::array<ctldl::Entry, 13> getRepeatingMtxEntriesA() {
  return {{
    {0,0},
    {4,0},
    {1,1},
    {4,1},
    {5,1},
    {2,2},
    {6,2},
    {3,3},
    {5,3},
    {6,3},
    {4,4},
    {5,5},
    {6,6},
  }};
}

constexpr std::array<ctldl::Entry, 6> getRepeatingMtxEntriesB() {
  return {{
    {4,0},
    {4,1},
    {5,1},
    {6,2},
    {5,3},
    {6,3},
  }};
}

constexpr std::array<std::size_t, getRepeatingMtxDim()>
getRepeatingMtxPermutation() {
  return {
    0,
    1,
    2,
    3,
    4,
    5,
    6,
  };
}
