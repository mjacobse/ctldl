#pragma once

#include <cstddef>
#include <type_traits>
#include <utility>

#include <ctldl/utility/make_index_sequence.hpp>

namespace ctldl {

template <typename F, std::size_t... Is>
void unroll(F&& f, std::index_sequence<Is...>) {
  (f(std::integral_constant<std::size_t, Is>{}), ...);
}

template <std::size_t i_begin, std::size_t i_end, typename F>
void unroll(F&& f) {
  unroll(std::forward<F>(f), makeIndexSequence<i_begin, i_end>());
}

template <std::size_t i_begin, std::size_t i_end, typename F>
void unrollReversed(F&& f) {
  unroll(std::forward<F>(f), makeIndexSequenceReversed<i_begin, i_end>());
}

}  // namespace ctldl
