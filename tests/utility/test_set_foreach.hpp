#pragma once

#include "tests/utility/test_set.hpp"

#include <functional>
#include <meta>

namespace ctldl {

template <template <class...> class Callable, std::meta::info test_set_types,
          typename... Args>
void foreach (const TestSetValues<Args...>& values_set) {
  template for (constexpr auto& types : [:test_set_types:]) {
    for (const auto& values : values_set) {
      constexpr std::meta::info callable_type =
          std::meta::substitute(^^Callable, types);
      using CallableType = [:callable_type:];
      std::apply(CallableType{}, values);
    }
  }
}

}  // namespace ctldl
