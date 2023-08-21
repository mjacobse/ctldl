#pragma once

#include <boost/test/unit_test.hpp>

#define CTLDL_TEST_STATIC_FIRST_ARG(X, ...) X
#define CTLDL_TEST_STATIC(...)                                           \
  {                                                                      \
    constexpr bool condition = CTLDL_TEST_STATIC_FIRST_ARG(__VA_ARGS__); \
    static_cast<void>(condition);                                        \
    BOOST_TEST(__VA_ARGS__);                                             \
  }
