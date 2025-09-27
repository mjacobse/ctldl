#pragma once

#include <source_location>

namespace ctldl {

void contractHandlePreconditionViolation(std::source_location source_location);
void contractHandlePostconditionViolation(std::source_location source_location);
void contractHandleAssertionViolation(std::source_location source_location);

inline void pre(bool condition, const std::source_location source_location =
                                    std::source_location{}) {
  if (condition) {
    return;
  }
  contractHandlePreconditionViolation(source_location);
}

inline void post(bool condition, const std::source_location source_location =
                                     std::source_location{}) {
  if (condition) {
    return;
  }
  contractHandlePostconditionViolation(source_location);
}

inline void contract_assert(
    bool condition,
    const std::source_location source_location = std::source_location{}) {
  if (condition) {
    return;
  }
  contractHandleAssertionViolation(source_location);
}

}  // namespace ctldl
