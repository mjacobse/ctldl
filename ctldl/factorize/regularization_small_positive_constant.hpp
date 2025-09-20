#pragma once

#include <cmath>
#include <cstddef>
#include <limits>

namespace ctldl {

struct RegularizationSmallPositiveConstant {
  template <class Value>
  Value regularize(const Value Di, [[maybe_unused]] std::size_t i) const {
    const auto pivot_tolerance = std::pow(std::numeric_limits<Value>::epsilon(),
                                          static_cast<Value>(0.625));
    const auto regularization_constant = std::pow(
        std::numeric_limits<Value>::epsilon(), static_cast<Value>(0.5));
    if (std::abs(Di) < pivot_tolerance) {
      return regularization_constant;
    }
    return Di;
  }
};

}  // namespace ctldl
