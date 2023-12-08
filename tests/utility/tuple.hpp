#pragma once

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t i, class T>
struct TupleElement {
  T t;
};

template <std::size_t i, class T>
const T& get(const TupleElement<i, T>& tuple) {
  return tuple.t;
}

template <class IndexSequence, class... Args>
struct TupleHelper {};

template <std::size_t... Is, class... Args>
struct TupleHelper<std::index_sequence<Is...>, Args...>
    : public TupleElement<Is, Args>... {};

/**
 * Custom simplified Tuple type that is much better on compilation times than
 * std::tuple.
 */
template <class... Args>
struct Tuple
    : public TupleHelper<std::make_index_sequence<sizeof...(Args)>, Args...> {};

template <class... L, class... R, std::size_t... IsL, std::size_t... IsR>
auto tupleCat(const Tuple<L...>& l, const Tuple<R...>& r,
              std::index_sequence<IsL...>, std::index_sequence<IsR...>) {
  return Tuple<L..., R...>{{{get<IsL>(l)}..., {get<IsR>(r)}...}};
}

template <class... L, class... R>
auto tupleCat(const Tuple<L...>& l, const Tuple<R...>& r) {
  return tupleCat(l, r, std::make_index_sequence<sizeof...(L)>(),
                  std::make_index_sequence<sizeof...(R)>());
}

template <class Callable, class TupleT, std::size_t... Is>
void apply(Callable&& callable, const TupleT& args,
           std::index_sequence<Is...>) {
  callable(get<Is>(args)...);
}

template <class Callable, class... Args>
void apply(Callable&& callable, const Tuple<Args...>& args) {
  apply(callable, args, std::make_index_sequence<sizeof...(Args)>());
}

}  // namespace ctldl
