#pragma once

#include <ctldl/sparsity/is_nonzero_info.hpp>

namespace ctldl {

template <class Sparsity>
constexpr auto getIsNonzeroInfo() {
  IsNonzeroInfo<Sparsity::num_rows, Sparsity::num_cols> is_nonzero;
  for (const auto entry : Sparsity::entries) {
    is_nonzero[entry.row_index][entry.col_index] = true;
  }
  return is_nonzero;
}

}  // namespace ctldl
