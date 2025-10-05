#pragma once

#include <ctldl/sparsity/sparsity.hpp>

namespace ctldl {

constexpr bool isLowerTriangle(const SparsityView sparsity) {
  for (const auto& entry : sparsity.entries()) {
    if (entry.row_index <= entry.col_index) {
      return false;
    }
  }
  return true;
}

}  // namespace ctldl
