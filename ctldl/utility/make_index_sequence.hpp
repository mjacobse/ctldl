#pragma once

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t Offset, std::size_t... Is>
constexpr std::index_sequence<(Offset + Is)...> addOffset(
    std::index_sequence<Is...>) {
  return {};
}

template <std::size_t Begin, std::size_t End>
constexpr auto makeIndexSequence() {
  return addOffset<Begin>(std::make_index_sequence<End - Begin>{});
}

template <std::size_t... Is>
constexpr std::index_sequence<(sizeof...(Is) - 1 - Is)...> reverse(
    std::index_sequence<Is...>) {
  return {};
}

template <std::size_t Begin, std::size_t End>
constexpr auto makeIndexSequenceReversed() {
  return addOffset<Begin>(reverse(std::make_index_sequence<End - Begin>{}));
}

}  // namespace ctldl
