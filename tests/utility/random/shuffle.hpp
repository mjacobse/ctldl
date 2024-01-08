#pragma once

#include "tests/utility/random/uniform_int_distribution.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>

namespace ctldl {

template <class RandomIt, class Generator>
constexpr RandomIt shuffle(const RandomIt first, const RandomIt last,
                           Generator& generator) {
  using DifferenceType =
      typename std::iterator_traits<RandomIt>::difference_type;
  using Distribution = UniformIntDistribution<std::ptrdiff_t>;

  DifferenceType num_remaining = last - first;
  for (auto it = first; it < last; ++it) {
    num_remaining -= 1;
    Distribution distribution(0, num_remaining);
    const DifferenceType i = distribution(generator);
    if (i != DifferenceType{0}) {
      std::iter_swap(it, it + i);
    }
  }
  return last;
}

}  // namespace ctldl
