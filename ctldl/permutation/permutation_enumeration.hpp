#pragma once

#include <ctldl/permutation/factorial.hpp>
#include <ctldl/permutation/permutation.hpp>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

namespace ctldl {

consteval std::vector<PermutationViewStructural> enumeratePermutationsStatic(
    const std::size_t dim) {
  auto current = std::views::iota(0uz, dim) | std::ranges::to<std::vector>();
  std::vector<PermutationViewStructural> ret;
  for (std::size_t i = 0; i < factorial(dim); ++i) {
    ret.push_back(PermutationDynamic(current));
    std::ranges::next_permutation(current);
  }
  return ret;
}

}  // namespace ctldl
