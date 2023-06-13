#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <span>
#include <vector>

namespace ctldl {

template <class Value, std::size_t block_dim>
std::vector<Value> flatten(
    const std::span<const std::array<Value, block_dim>> values_blocked) {
  std::vector<Value> values_flattened;
  for (const auto& block : values_blocked) {
    std::copy(block.cbegin(), block.cend(),
              std::back_inserter(values_flattened));
  }
  return values_flattened;
}

}  // namespace ctldl
