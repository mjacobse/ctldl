#pragma once

#include <cstddef>
#include <meta>
#include <ranges>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace ctldl {

// Structural version of std::span<const std::meta::info> so that we can use
// std::define_static_array for a range of such elements
struct SpanConstMetaInfo {
  consteval explicit(false)
      SpanConstMetaInfo(const std::span<const std::meta::info> init)
      : m_ptr(init.data()), m_size(init.size()) {}

  consteval auto begin() const { return m_ptr; }
  consteval auto end() const { return m_ptr + m_size; }

  const std::meta::info* m_ptr;
  std::size_t m_size;
};

class TestSetTypes {
 public:
  consteval explicit TestSetTypes(const std::ranges::range auto& init)
      : m_instances(std::define_static_array(init)) {}

  consteval auto begin() const { return m_instances.begin(); }
  consteval auto end() const { return m_instances.end(); }

 private:
  std::span<const SpanConstMetaInfo> m_instances;
};

template <class... Types>
consteval auto makeTypeArgument() {
  return TestSetTypes{std::views::transform(
      std::vector{^^Types...}, [](const std::meta::info& type) {
        return SpanConstMetaInfo{std::define_static_array(std::vector{type})};
      })};
}

template <class... Types>
consteval auto makeTypeArgumentFromTuple(std::tuple<Types...>) {
  return makeTypeArgument<Types...>();
}

consteval auto cartesianProductTestSets(const TestSetTypes& lhs,
                                        const TestSetTypes& rhs) {
  return TestSetTypes{std::views::cartesian_product(lhs, rhs) |
                      std::views::transform([](const auto& tuple) {
                        const auto& [... parts] = tuple;
                        return SpanConstMetaInfo{std::define_static_array(
                            std::views::concat(parts...))};
                      })};
}

consteval auto joinTestSets(const TestSetTypes& lhs, const TestSetTypes& rhs) {
  return TestSetTypes{std::views::concat(lhs, rhs)};
}

consteval auto zipTestSets(const TestSetTypes& lhs, const TestSetTypes& rhs) {
  return TestSetTypes{
      std::views::zip(lhs, rhs) | std::views::transform([](const auto& tuple) {
        const auto& [... parts] = tuple;
        return SpanConstMetaInfo{
            std::define_static_array(std::views::concat(parts...))};
      })};
}

consteval auto operator*(const TestSetTypes& lhs, const TestSetTypes& rhs) {
  return cartesianProductTestSets(lhs, rhs);
}

consteval auto operator+(const TestSetTypes& lhs, const TestSetTypes& rhs) {
  return joinTestSets(lhs, rhs);
}

consteval auto operator^(const TestSetTypes& lhs, const TestSetTypes& rhs) {
  return zipTestSets(lhs, rhs);
}

template <typename... Args>
class TestSetValues {
 public:
  TestSetValues(const std::ranges::range auto& init)
      : m_instances(init | std::ranges::to<std::vector>()) {}

  auto begin() const { return m_instances.begin(); }
  auto end() const { return m_instances.end(); }

 private:
  std::vector<std::tuple<Args...>> m_instances;
};

template <class T>
auto makeValueArgument(const std::initializer_list<T>& ts) {
  return TestSetValues<T>{
      std::views::transform(ts, [](const T& t) { return std::tuple<T>(t); })};
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto cartesianProductTestSets(const TestSetValues<ArgsLhs...>& lhs,
                              const TestSetValues<ArgsRhs...>& rhs) {
  return TestSetValues<ArgsLhs..., ArgsRhs...>{
      std::views::cartesian_product(lhs, rhs) |
      std::views::transform([](const auto& tuple) {
        const auto& [... parts] = tuple;
        return std::tuple_cat(parts...);
      })};
}

template <typename... Args>
auto joinTestSets(const TestSetValues<Args...>& lhs,
                  const TestSetValues<Args...>& rhs) {
  return TestSetValues{std::views::concat(lhs, rhs)};
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto zipTestSets(const TestSetValues<ArgsLhs...>& lhs,
                 const TestSetValues<ArgsRhs...>& rhs) {
  return TestSetValues{std::views::zip(lhs, rhs) |
                       std::views::transform([](const auto& tuple) {
                         const auto& [... parts] = tuple;
                         return std::tuple_cat(parts...);
                       })};
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto operator*(const TestSetValues<ArgsLhs...>& lhs,
               const TestSetValues<ArgsRhs...>& rhs) {
  return cartesianProductTestSets(lhs, rhs);
}

template <typename... Args>
auto operator+(const TestSetValues<Args...>& lhs,
               const TestSetValues<Args...>& rhs) {
  return joinTestSets(lhs, rhs);
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto operator^(const TestSetValues<ArgsLhs...>& lhs,
               const TestSetValues<ArgsRhs...>& rhs) {
  return zipTestSets(lhs, rhs);
}

}  // namespace ctldl
