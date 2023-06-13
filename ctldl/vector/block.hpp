#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace ctldl {

template <std::size_t block_dim, class Value>
std::vector<std::array<Value, block_dim>> block(
    const std::span<const Value> values_flat) {
  const auto num_blocks = std::size_t{values_flat.size() / block_dim};
  std::vector<std::array<Value, block_dim>> values_blocked(num_blocks);
  auto it_values_flat = std::cbegin(values_flat);
  for (auto& block : values_blocked) {
    for (auto& value : block) {
      value = *it_values_flat;
      ++it_values_flat;
    }
  }
  return values_blocked;
}

}  // namespace ctldl
