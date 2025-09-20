#pragma once

#include <concepts>
#include <cstddef>

namespace ctldl {

template <class T>
concept Regularization = requires(const T& t) {
  { t.regularize(1.0, std::size_t{0}) } -> std::convertible_to<double>;
};

}  // namespace ctldl
