#pragma once

#include "tests/utility/test_set.hpp"
#include "tests/utility/tuple.hpp"

#include <numeric>

namespace ctldl {

template <template <class...> class Callable, class Types>
struct Apply {};

template <template <class...> class Callable, class... Types>
struct Apply<Callable, Tuple<Types...>> {
  template <class Values>
  void operator()(const Values& values) {
    apply(Callable<Types...>{}, values);
  }
};

template <template <class...> class Callable, TestSet Set>
void foreach (const Set& set) {
  template for (constexpr auto i_types :
                std::views::iota(0uz, Set::sizeTypes())) {
    const auto size_values = set.template sizeValues<i_types>();
    for (std::size_t i_values = 0; i_values < size_values; ++i_values) {
      const auto data_point = set.template get<i_types>(i_values);
      Apply<Callable, typename decltype(data_point)::TestTypes>{}(
          data_point.test_values);
    }
  }
}

}  // namespace ctldl
