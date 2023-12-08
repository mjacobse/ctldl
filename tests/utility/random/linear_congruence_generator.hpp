#pragma once

#include <cstddef>
#include <limits>

namespace ctldl {

struct LinearCongruenceGenerator {
  std::size_t state;

  constexpr explicit LinearCongruenceGenerator(const std::size_t seed)
      : state(seed) {}

  static constexpr auto min() {
    return std::numeric_limits<decltype(state)>::min();
  }

  static constexpr auto max() {
    return std::numeric_limits<decltype(state)>::max();
  }

  constexpr std::size_t operator()() {
    state = state * 6364136223846793005 + 1;
    return state;
  }
};

}  // namespace ctldl
