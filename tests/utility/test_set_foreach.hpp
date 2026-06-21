#pragma once

#include "tests/utility/test_set.hpp"

#include <functional>
#include <meta>

namespace ctldl {

template <std::meta::info CallableTemplate, std::meta::info test_set_template,
          typename... Args>
void foreach (const TestSetFunction<Args...>& test_set_function) {
  template for (constexpr auto& args_template : [:test_set_template:]) {
    for (const auto& args_function : test_set_function) {
      constexpr std::meta::info callable =
          std::meta::substitute(CallableTemplate, args_template);
      using Callable = [:callable:];
      std::apply(Callable{}, args_function);
    }
  }
}

}  // namespace ctldl
