#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

constexpr int getRepeatingMtxDim() { return 16; }

constexpr std::array<ctldl::Entry, 28> getRepeatingMtxEntriesA() {
  return {{{0, 0},   {1, 0},   {1, 1},   {2, 2},   {3, 2},   {3, 3},
           {4, 4},   {5, 4},   {5, 5},   {6, 4},   {6, 5},   {6, 6},
           {7, 4},   {7, 5},   {7, 6},   {7, 7},   {8, 8},   {9, 8},
           {9, 9},   {10, 10}, {11, 10}, {11, 11}, {12, 12}, {13, 12},
           {13, 13}, {14, 14}, {15, 14}, {15, 15}}};
};

constexpr std::array<ctldl::Entry, 26> getRepeatingMtxEntriesB() {
  return {{{0, 1},   {1, 1},   {2, 2},   {3, 2},   {2, 3},   {3, 3},   {4, 5},
           {5, 5},   {6, 5},   {7, 5},   {4, 7},   {5, 7},   {6, 7},   {7, 7},
           {8, 9},   {9, 9},   {10, 11}, {11, 11}, {12, 12}, {13, 12}, {12, 13},
           {13, 13}, {14, 14}, {15, 14}, {14, 15}, {15, 15}}};
}

constexpr std::array<std::size_t, getRepeatingMtxDim()>
getRepeatingMtxPermutation() {
  return {0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15};
}
