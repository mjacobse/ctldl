#pragma once

#include <cstddef>

namespace ctldl {

constexpr std::size_t factorial(const std::size_t n) {
  std::size_t value = 1;
  for (std::size_t i = value + 1; i <= n; ++i) {
    value *= i;
  }
  return value;
}

}  // namespace ctldl
