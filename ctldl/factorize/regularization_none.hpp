#pragma once

#include <cstddef>

namespace ctldl {

struct RegularizationNone {
  template <class Value>
  constexpr Value regularize(const Value Di,
                             [[maybe_unused]] std::size_t i) const {
    return Di;
  }
};

}  // namespace ctldl
