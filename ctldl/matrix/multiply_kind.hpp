#pragma once

#include <type_traits>

namespace ctldl {

enum class MultiplyKind {
  Normal = 1,
  Transposed = 2,
  Symmetric = Normal | Transposed,
};

constexpr bool operator&(const MultiplyKind lhs, const MultiplyKind rhs) {
  using T = std::underlying_type_t<MultiplyKind>;
  return static_cast<T>(lhs) & static_cast<T>(rhs);
}

}  // namespace ctldl
