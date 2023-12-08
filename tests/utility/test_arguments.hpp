#pragma once

#include "tests/utility/tuple.hpp"

#include <utility>

namespace ctldl {

template <class TypeTuple, class ValueTuple>
struct TestArguments {
  using TestTypes = TypeTuple;
  ValueTuple test_values;
};

template <class TypeTupleLhs, class ValueTupleLhs, class TypeTupleRhs,
          class ValueTupleRhs>
constexpr auto catTestArgs(
    const TestArguments<TypeTupleLhs, ValueTupleLhs>& lhs,
    const TestArguments<TypeTupleRhs, ValueTupleRhs>& rhs) {
  using TypeTupleCat = decltype(tupleCat(std::declval<TypeTupleLhs>(),
                                         std::declval<TypeTupleRhs>()));
  using ValueTupleCat = decltype(tupleCat(std::declval<ValueTupleLhs>(),
                                          std::declval<ValueTupleRhs>()));
  return TestArguments<TypeTupleCat, ValueTupleCat>{
      tupleCat(lhs.test_values, rhs.test_values)};
}

}  // namespace ctldl
