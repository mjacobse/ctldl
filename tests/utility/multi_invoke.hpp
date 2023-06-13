#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

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

template <class Types, class NewType>
struct TupleAdd {};

template <class... Types, class NewType>
struct TupleAdd<std::tuple<Types...>, NewType> {
  using type = std::tuple<Types..., NewType>;
};

template <class Types, class NewType>
using TupleAddT = typename TupleAdd<Types, NewType>::type;

template <template <class...> class Callable, std::size_t i_type,
          class TypeOptionsTuple, class ChosenTypeTuple, std::size_t... Is,
          class ChosenValueTuple>
void multiInvokeIterateTypeOptions(const ChosenValueTuple& chosen_value_tuple,
                                   std::index_sequence<Is...>) {
  using TypeOptionsTupleExtended = TupleAddT<TypeOptionsTuple, std::tuple<>>;
  if constexpr (i_type < std::tuple_size_v<TypeOptionsTuple>) {
    using TypeOptions = std::tuple_element_t<i_type, TypeOptionsTuple>;
    using TypeOptionsNext =
        std::tuple_element_t<i_type + 1, TypeOptionsTupleExtended>;
    constexpr auto num_type_options_next = std::tuple_size_v<TypeOptionsNext>;
    (multiInvokeIterateTypeOptions<
         Callable, i_type + 1, TypeOptionsTuple,
         TupleAddT<ChosenTypeTuple, std::tuple_element_t<Is, TypeOptions>>>(
         chosen_value_tuple, std::make_index_sequence<num_type_options_next>()),
     ...);
  } else {
    Apply<Callable, ChosenTypeTuple>{}(chosen_value_tuple);
  }
}

template <template <class...> class Callable, std::size_t i_value,
          class TypeOptionsTuple, class ChosenValueTuple,
          class ValueOptionsTuple>
void multiInvokeIterateValueOptions(
    const ValueOptionsTuple& value_options_tuple,
    const ChosenValueTuple& chosen_value_tuple) {
  if constexpr (i_value < std::tuple_size_v<ValueOptionsTuple>) {
    const auto& value_options = std::get<i_value>(value_options_tuple);
    for (const auto& value : value_options) {
      multiInvokeIterateValueOptions<Callable, i_value + 1, TypeOptionsTuple>(
          value_options_tuple,
          std::tuple_cat(chosen_value_tuple, std::make_tuple(value)));
    }
  } else {
    if constexpr (std::tuple_size_v<TypeOptionsTuple> == 0) {
      std::apply(Callable<>{}, chosen_value_tuple);
    } else {
      using TypeOptions = std::tuple_element_t<0, TypeOptionsTuple>;
      constexpr auto num_type_options = std::tuple_size_v<TypeOptions>;
      multiInvokeIterateTypeOptions<Callable, 0, TypeOptionsTuple,
                                    std::tuple<>>(
          chosen_value_tuple, std::make_index_sequence<num_type_options>());
    }
  }
}

template <template <class...> class Callable, class TypeOptionsTuple,
          class ValueOptionsTuple>
void multiInvoke(const ValueOptionsTuple& value_options_tuple) {
  multiInvokeIterateValueOptions<Callable, 0, TypeOptionsTuple>(
      value_options_tuple, std::tuple<>{});
}

}  // namespace ctldl
