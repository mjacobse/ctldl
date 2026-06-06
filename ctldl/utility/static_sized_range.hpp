#pragma once

#include <cstddef>
#include <ranges>
#include <type_traits>
#include <utility>

namespace ctldl {

template <class T>
concept static_sized_range = std::ranges::sized_range<T> && requires(T& t) {
  std::integral_constant<std::size_t, std::ranges::size(t)>{};
};

template <static_sized_range T>
constexpr auto range_static_size_v = decltype([](T& t) {
  return std::integral_constant<std::size_t, std::ranges::size(t)>{};
}(std::declval<T&>()))::value;

}  // namespace ctldl
