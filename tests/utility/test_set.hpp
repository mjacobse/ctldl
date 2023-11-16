#pragma once

#include "tests/utility/is_specialization_of.hpp"
#include "tests/utility/test_arguments.hpp"

#include <algorithm>
#include <cassert>
#include <concepts>
#include <memory>
#include <tuple>
#include <type_traits>

namespace ctldl {

template <typename T>
concept TestSet = requires(T t) {
  // static function sizeTypes returns the number of different type argument
  // combinations contained in the test set
  { T::sizeTypes() } -> std::convertible_to<std::size_t>;
  // member function sizeValues<i> returns the number of different value
  // argument combinations contained in the test set for the i-th type argument
  // combination
  { t.template sizeValues<0>() } -> std::convertible_to<std::size_t>;
  // member function get<i>(j) returns the i-th type argument combination
  // together with the j-th value argument combination
  requires is_specialization_of_v<decltype(t.template get<0>(0)),
                                  TestArguments>;
};

template <TestSet Left, TestSet Right>
class TestSetCombinationCartesian {
 public:
  TestSetCombinationCartesian(const Left& left, const Right& right)
      : m_left(left), m_right(right) {}

  static constexpr auto sizeTypes() {
    return Left::sizeTypes() * Right::sizeTypes();
  }

  template <std::size_t i_types>
  constexpr auto sizeValues() const {
    constexpr auto i_left = indexTypesLeft<i_types>();
    constexpr auto i_right = indexTypesRight<i_types>();
    return m_left.template sizeValues<i_left>() *
           m_right.template sizeValues<i_right>();
  }

  template <std::size_t i_types>
  constexpr auto get(const std::size_t i_values) const {
    constexpr auto i_types_left = indexTypesLeft<i_types>();
    constexpr auto i_types_right = indexTypesRight<i_types>();
    const auto size_values = m_right.template sizeValues<i_types_right>();
    const auto i_values_left = i_values / size_values;
    const auto i_values_right = i_values % size_values;
    return catTestArgs(m_left.template get<i_types_left>(i_values_left),
                          m_right.template get<i_types_right>(i_values_right));
  }

 private:
  Left m_left;
  Right m_right;

  template <std::size_t i_types>
  static constexpr auto indexTypesLeft() {
    return i_types / Right::sizeTypes();
  }

  template <std::size_t i_types>
  static constexpr auto indexTypesRight() {
    return i_types % Right::sizeTypes();
  }
};

template <TestSet Left, TestSet Right>
class TestSetCombinationJoin {
 public:
  TestSetCombinationJoin(const Left& left, const Right& right)
      : m_left(left), m_right(right) {}

  static constexpr auto sizeTypes() {
    return Left::sizeTypes() + Right::sizeTypes();
  }

  template <std::size_t i_types>
  constexpr auto sizeValues() const {
    if constexpr (i_types < Left::sizeTypes()) {
      constexpr auto i_types_left = indexTypesLeft<i_types>();
      return m_left.template sizeValues<i_types_left>();
    } else {
      constexpr auto i_types_right = indexTypesRight<i_types>();
      return m_right.template sizeValues<i_types_right>();
    }
  }

  template <std::size_t i_types>
  constexpr auto get(const std::size_t i_values) const {
    if constexpr (i_types < Left::sizeTypes()) {
      constexpr auto i_types_left = indexTypesLeft<i_types>();
      return m_left.template get<i_types_left>(i_values);
    } else {
      constexpr auto i_types_right = indexTypesRight<i_types>();
      return m_right.template get<i_types_right>(i_values);
    }
  }

 private:
  Left m_left;
  Right m_right;

  template <std::size_t i_types>
  static constexpr auto indexTypesLeft() {
    return i_types;
  }

  template <std::size_t i_types>
  static constexpr auto indexTypesRight() {
    return i_types - Left::sizeTypes();
  }
};

template <TestSet Left, TestSet Right>
class TestSetCombinationZip {
 public:
  TestSetCombinationZip(const Left& left, const Right& right)
      : m_left(left), m_right(right) {}

  static_assert(Left::sizeTypes() == Right::sizeTypes());
  static constexpr auto sizeTypes() {
    return Left::sizeTypes();
  }

  template <std::size_t i_types>
  constexpr auto sizeValues() const {
    assert(m_left.template sizeValues<i_types>() ==
           m_right.template sizeValues<i_types>());
    return m_left.template sizeValues<i_types>();
  }

  template <std::size_t i_types>
  constexpr auto get(const std::size_t i_values) const {
    return catTestArgs(m_left.template get<i_types>(i_values),
                          m_right.template get<i_types>(i_values));
  }

 private:
  Left m_left;
  Right m_right;
};

template <class Types>
struct TypeArgument {
  static constexpr auto sizeTypes() {
    return std::tuple_size_v<Types>;
  }

  template <std::size_t i_types>
  static constexpr auto sizeValues() {
    return std::size_t{1};
  }

  template <std::size_t i_types>
  static constexpr auto get(const std::size_t i_values) {
    assert(i_values == 0);
    static_cast<void>(i_values);
    using TestTypes = std::tuple<std::tuple_element_t<i_types, Types>>;
    using TestValues = std::tuple<>;
    return TestArguments<TestTypes, TestValues>{};
  }
};

template <class... Types>
auto makeTypeArgument() {
  return TypeArgument<std::tuple<Types...>>{};
}

template <class T>
class ValueArgument {
 public:
  ValueArgument(const std::initializer_list<T> init_list)
      : m_values([init_list] {
          auto values = std::make_unique<T[]>(init_list.size());
          std::copy(std::cbegin(init_list), std::cend(init_list), values.get());
          return values;
        }()),
        m_num_values(init_list.size()) {}

  static constexpr auto sizeTypes() { return 1; }

  template <std::size_t i_types>
  constexpr auto sizeValues() const {
    static_assert(i_types == 0);
    return m_num_values;
  }

  template <std::size_t i_types>
  constexpr auto get(const std::size_t i_values) const {
    static_assert(i_types == 0);
    assert(i_values < m_num_values);
    return TestArguments<std::tuple<>, std::tuple<T>>{
        std::make_tuple(m_values[static_cast<std::ptrdiff_t>(i_values)])};
  }

 private:
  std::shared_ptr<const T[]> m_values;
  std::size_t m_num_values;
};

template <class T>
auto makeValueArgument(const std::initializer_list<T> ts) {
  return ValueArgument<T>{ts};
}

template <TestSet Lhs, TestSet Rhs>
auto cartesianProductTestSets(const Lhs& lhs, const Rhs& rhs) {
  return TestSetCombinationCartesian<Lhs, Rhs>{lhs, rhs};
}

template <TestSet Lhs, TestSet Rhs>
auto joinTestSets(const Lhs& lhs, const Rhs& rhs) {
  return TestSetCombinationJoin<Lhs, Rhs>{lhs, rhs};
}

template <TestSet Lhs, TestSet Rhs>
auto zipTestSets(const Lhs& lhs, const Rhs& rhs) {
  return TestSetCombinationZip<Lhs, Rhs>{lhs, rhs};
}

template <TestSet Lhs, TestSet Rhs>
auto operator*(const Lhs& lhs, const Rhs& rhs) {
  return cartesianProductTestSets(lhs, rhs);
}

template <TestSet Lhs, TestSet Rhs>
auto operator+(const Lhs& lhs, const Rhs& rhs) {
  return joinTestSets(lhs, rhs);
}

template <TestSet Lhs, TestSet Rhs>
auto operator^(const Lhs& lhs, const Rhs& rhs) {
  return zipTestSets(lhs, rhs);
}

}  // namespace ctldl
