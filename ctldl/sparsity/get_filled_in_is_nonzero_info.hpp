#pragma once

#include <ctldl/sparsity/is_nonzero_info.hpp>

namespace ctldl {

template <int dim>
constexpr auto getFilledInIsNonzeroInfo(
    const IsNonzeroInfo<dim, dim> is_nonzero_original_matrix) {
  IsNonzeroInfo<dim, dim> is_nonzero;
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < i; ++j) {
      is_nonzero[i][j] = is_nonzero_original_matrix[i][j];
      for (int k = 0; k < j; ++k) {
        if (is_nonzero[i][k] && is_nonzero[j][k]) {
          is_nonzero[i][j] = true;
        }
      }
    }
  }
  return is_nonzero;
}

}  // namespace ctldl
