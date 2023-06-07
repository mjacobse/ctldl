#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

constexpr int getRepeatingMtxDim() { return 20; }

constexpr std::array<ctldl::Entry, 45> getRepeatingMtxEntriesA() {
  return {{
    {0,0},
    {11,0},
    {1,1},
    {11,1},
    {12,1},
    {2,2},
    {5,2},
    {7,2},
    {13,2},
    {14,2},
    {3,3},
    {13,3},
    {14,3},
    {4,4},
    {12,4},
    {14,4},
    {15,4},
    {5,5},
    {7,5},
    {14,5},
    {16,5},
    {6,6},
    {16,6},
    {17,6},
    {7,7},
    {14,7},
    {17,7},
    {18,7},
    {8,8},
    {19,8},
    {9,9},
    {15,9},
    {19,9},
    {10,10},
    {18,10},
    {19,10},
    {11,11},
    {12,12},
    {13,13},
    {14,14},
    {15,15},
    {16,16},
    {17,17},
    {18,18},
    {19,19},
  }};
}

constexpr std::array<ctldl::Entry, 22> getRepeatingMtxEntriesB() {
  return {{
    {11,0},
    {11,1},
    {12,1},
    {13,2},
    {14,2},
    {13,3},
    {14,3},
    {12,4},
    {14,4},
    {15,4},
    {14,5},
    {16,5},
    {16,6},
    {17,6},
    {14,7},
    {17,7},
    {18,7},
    {19,8},
    {15,9},
    {19,9},
    {18,10},
    {19,10},
  }};
}

constexpr std::array<std::size_t, getRepeatingMtxDim()>
getRepeatingMtxPermutation() {
  return {
    6,
    1,
    0,
    4,
    5,
    10,
    2,
    9,
    3,
    11,
    8,
    12,
    15,
    18,
    19,
    14,
    13,
    7,
    17,
    16,
  };
}
