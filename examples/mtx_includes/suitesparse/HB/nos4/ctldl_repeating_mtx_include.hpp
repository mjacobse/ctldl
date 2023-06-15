#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

constexpr int getRepeatingMtxDim() { return 10; }

constexpr std::array<ctldl::Entry, 16> getRepeatingMtxEntriesA() {
  return {{
    {0,0},
    {1,0},
    {2,0},
    {1,1},
    {2,2},
    {4,2},
    {3,3},
    {4,4},
    {6,4},
    {5,5},
    {6,6},
    {8,6},
    {7,7},
    {8,8},
    {9,8},
    {9,9},
  }};
}

constexpr std::array<ctldl::Entry, 21> getRepeatingMtxEntriesB() {
  return {{
    {2,0},
    {3,0},
    {1,1},
    {2,1},
    {3,1},
    {3,3},
    {2,4},
    {3,4},
    {6,4},
    {7,4},
    {2,5},
    {3,5},
    {5,5},
    {6,5},
    {7,5},
    {7,7},
    {6,8},
    {7,8},
    {6,9},
    {7,9},
    {9,9},
  }};
}

constexpr std::array<std::size_t, getRepeatingMtxDim()>
getRepeatingMtxPermutation() {
  return {7, 8, 0, 4, 3, 2, 6, 5, 9, 1};
}
