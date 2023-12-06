#pragma once

#include <array>
#include <cstddef>

namespace ctldl {

// This is to workaround https://github.com/llvm/llvm-project/issues/74375
// We could just use {} on the array directly, but then we lose the benefit of
// diagnostics in case the actual initialization that we want to do misses any
// entries.
template <class T, std::size_t N>
constexpr void fixInitIfZeroLengthArray(std::array<T, N>& arr) {
  if constexpr (N == 0) {
    arr = std::array<T, N>{};
  }
}

}  // namespace ctldl
