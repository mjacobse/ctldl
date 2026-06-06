#pragma once

#include <array>
#include <cstddef>
#include <ranges>
#include <type_traits>

namespace ctldl {

template <class T>
struct range_static_size;

template <class T, std::size_t n>
struct range_static_size<std::array<T, n>>
    : std::integral_constant<std::size_t, n> {};

template <class T, std::size_t n>
struct range_static_size<T[n]> : std::integral_constant<std::size_t, n> {};

template <class T>
inline constexpr std::size_t range_static_size_v = range_static_size<T>::value;

template <class T>
concept static_sized_range = std::ranges::sized_range<T> && requires {
  std::integral_constant<std::size_t, range_static_size_v<T>>{};
};

}  // namespace ctldl
