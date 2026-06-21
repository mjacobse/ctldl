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

class TestSetTemplate {
 public:
  consteval explicit TestSetTemplate(const std::ranges::range auto& init)
      : m_instances(std::define_static_array(init)) {}

  consteval auto begin() const { return m_instances.begin(); }
  consteval auto end() const { return m_instances.end(); }

 private:
  std::span<const SpanConstMetaInfo> m_instances;
};

template <class... Types>
consteval auto makeTemplateArgument() {
  return TestSetTemplate{std::views::transform(
      std::vector{^^Types...}, [](const std::meta::info& type) {
        return SpanConstMetaInfo{std::define_static_array(std::vector{type})};
      })};
}

consteval auto makeTemplateArgument(const std::ranges::range auto& r) {
  return TestSetTemplate{std::views::transform(r, [](const auto& value) {
    return SpanConstMetaInfo{std::define_static_array(
        std::vector{std::meta::reflect_constant(value)})};
  })};
}

consteval auto cartesianProductTestSets(const TestSetTemplate& lhs,
                                        const TestSetTemplate& rhs) {
  return TestSetTemplate{std::views::cartesian_product(lhs, rhs) |
                         std::views::transform([](const auto& tuple) {
                           const auto& [... parts] = tuple;
                           return SpanConstMetaInfo{std::define_static_array(
                               std::views::concat(parts...))};
                         })};
}

consteval auto joinTestSets(const TestSetTemplate& lhs,
                            const TestSetTemplate& rhs) {
  return TestSetTemplate{std::views::concat(lhs, rhs)};
}

consteval auto zipTestSets(const TestSetTemplate& lhs,
                           const TestSetTemplate& rhs) {
  return TestSetTemplate{
      std::views::zip(lhs, rhs) | std::views::transform([](const auto& tuple) {
        const auto& [... parts] = tuple;
        return SpanConstMetaInfo{
            std::define_static_array(std::views::concat(parts...))};
      })};
}

consteval auto operator*(const TestSetTemplate& lhs,
                         const TestSetTemplate& rhs) {
  return cartesianProductTestSets(lhs, rhs);
}

consteval auto operator+(const TestSetTemplate& lhs,
                         const TestSetTemplate& rhs) {
  return joinTestSets(lhs, rhs);
}

consteval auto operator^(const TestSetTemplate& lhs,
                         const TestSetTemplate& rhs) {
  return zipTestSets(lhs, rhs);
}

template <typename... Args>
class TestSetFunction {
 public:
  TestSetFunction(const std::ranges::range auto& init)
      : m_instances(init | std::ranges::to<std::vector>()) {}

  auto begin() const { return m_instances.begin(); }
  auto end() const { return m_instances.end(); }

 private:
  std::vector<std::tuple<Args...>> m_instances;
};

template <class T>
auto makeFunctionArgument(const std::initializer_list<T>& ts) {
  return TestSetFunction<T>{
      std::views::transform(ts, [](const T& t) { return std::tuple<T>(t); })};
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto cartesianProductTestSets(const TestSetFunction<ArgsLhs...>& lhs,
                              const TestSetFunction<ArgsRhs...>& rhs) {
  return TestSetFunction<ArgsLhs..., ArgsRhs...>{
      std::views::cartesian_product(lhs, rhs) |
      std::views::transform([](const auto& tuple) {
        const auto& [... parts] = tuple;
        return std::tuple_cat(parts...);
      })};
}

template <typename... Args>
auto joinTestSets(const TestSetFunction<Args...>& lhs,
                  const TestSetFunction<Args...>& rhs) {
  return TestSetFunction{std::views::concat(lhs, rhs)};
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto zipTestSets(const TestSetFunction<ArgsLhs...>& lhs,
                 const TestSetFunction<ArgsRhs...>& rhs) {
  return TestSetFunction{std::views::zip(lhs, rhs) |
                         std::views::transform([](const auto& tuple) {
                           const auto& [... parts] = tuple;
                           return std::tuple_cat(parts...);
                         })};
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto operator*(const TestSetFunction<ArgsLhs...>& lhs,
               const TestSetFunction<ArgsRhs...>& rhs) {
  return cartesianProductTestSets(lhs, rhs);
}

template <typename... Args>
auto operator+(const TestSetFunction<Args...>& lhs,
               const TestSetFunction<Args...>& rhs) {
  return joinTestSets(lhs, rhs);
}

template <typename... ArgsLhs, typename... ArgsRhs>
auto operator^(const TestSetFunction<ArgsLhs...>& lhs,
               const TestSetFunction<ArgsRhs...>& rhs) {
  return zipTestSets(lhs, rhs);
}

}  // namespace ctldl
