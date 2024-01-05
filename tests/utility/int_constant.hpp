#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

namespace ctldl {

template <int value_>
struct IntConstant {
  static constexpr auto value = value_;
};

template <std::size_t... Is>
constexpr auto makeIntConstantSequence(std::index_sequence<Is...>) {
  return std::tuple<IntConstant<int{Is}>...>{};
}

template <std::size_t num>
constexpr auto makeIntConstantSequence() {
  return makeIntConstantSequence(std::make_index_sequence<num>());
}

}  // namespace ctldl
