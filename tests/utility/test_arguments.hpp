#pragma once

#include <tuple>
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
  using TypeTupleCat = decltype(std::tuple_cat(std::declval<TypeTupleLhs>(),
                                               std::declval<TypeTupleRhs>()));
  using ValueTupleCat = decltype(std::tuple_cat(std::declval<ValueTupleLhs>(),
                                                std::declval<ValueTupleRhs>()));
  return TestArguments<TypeTupleCat, ValueTupleCat>{
      std::tuple_cat(lhs.test_values, rhs.test_values)};
}

}  // namespace ctldl
