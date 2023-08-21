#pragma once

namespace ctldl {

template <class Sparsity>
constexpr bool isLowerTriangle(const Sparsity& sparsity) {
  for (const auto& entry : sparsity.entries) {
    if (entry.row_index <= entry.col_index) {
      return false;
    }
  }
  return true;
}

}  // namespace ctldl
