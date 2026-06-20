#pragma once

#include <algorithm>
#include <initializer_list>
#include <ranges>

namespace ctldl {

template <typename T>
constexpr bool all_equal(const std::initializer_list<T> values) {
  return std::ranges::adjacent_find(values, std::ranges::not_equal_to{}) ==
         values.end();
}

}  // namespace ctldl
