#pragma once

#include "tests/utility/test_set.hpp"

#include <numeric>
#include <tuple>

namespace ctldl {

template <template <class...> class Callable, class Types>
struct Apply {};

template <template <class...> class Callable, class... Types>
struct Apply<Callable, std::tuple<Types...>> {
  template <class Values>
  void operator()(const Values& values) {
    std::apply(Callable<Types...>{}, values);
  }
};

template <template <class...> class Callable, std::size_t i_types, TestSet Set>
void foreachFold(const Set& set) {
  const auto size_values = set.template sizeValues<i_types>();
  for (std::size_t i_values = 0; i_values < size_values; ++i_values) {
    const auto data_point = set.template get<i_types>(i_values);
    Apply<Callable, typename decltype(data_point)::TestTypes>{}(
        data_point.test_values);
  }
}

template <template <class...> class Callable, std::size_t offset,
          std::size_t... Is, TestSet Set>
void foreachFold(const Set& set, std::index_sequence<Is...>) {
  (foreachFold<Callable, offset + Is>(set), ...);
}

template <template <class...> class Callable, std::size_t i_begin,
          std::size_t i_end, TestSet Set>
void foreachIter(const Set& set) {
  // Split the set of type argument combinations in half until small enough,
  // then do a fold expression. This avoids hitting compiler nesting limit for
  // fold expressions and also avoids deep recursive template instantiation.
  constexpr auto num = i_end - i_begin;
  constexpr auto limit_fold = std::size_t{8};
  if constexpr (num <= limit_fold) {
    foreachFold<Callable, i_begin>(set, std::make_index_sequence<num>());
  } else {
    constexpr auto i_mid = std::midpoint(i_begin, i_end);
    foreachIter<Callable, i_begin, i_mid>(set);
    foreachIter<Callable, i_mid, i_end>(set);
  }
}

template <template <class...> class Callable, TestSet Set>
void foreach (const Set& set) {
  foreachIter<Callable, 0, Set::sizeTypes()>(set);
}

}  // namespace ctldl
