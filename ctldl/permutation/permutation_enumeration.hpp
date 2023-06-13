#pragma once

#include <ctldl/permutation/factorial.hpp>
#include <ctldl/permutation/permutation.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <utility>

namespace ctldl {

template <std::size_t dim, std::size_t i>
struct EnumeratePermutationsHelper {
  static constexpr std::array<std::size_t, dim> permutation = [] {
    std::array<std::size_t, dim> ret;
    if constexpr (i == 0) {
      std::iota(ret.begin(), ret.end(), 0);
    } else {
      ret = EnumeratePermutationsHelper<dim, i - 1>::permutation;
      std::next_permutation(ret.begin(), ret.end());
    }
    return ret;
  }();
};

template <std::size_t dim, std::size_t... Is>
constexpr auto enumeratePermutations(std::index_sequence<Is...>) {
  return std::make_tuple(EnumeratePermutationsHelper<dim, Is>{}...);
}

template <std::size_t dim>
using PermutationEnumeration = decltype(enumeratePermutations<dim>(
    std::make_index_sequence<factorial(dim)>()));

}  // namespace ctldl
