#pragma once

#include <ctldl/sparsity/is_nonzero_info.hpp>

#include <cstddef>

namespace ctldl {

template <std::size_t dim>
constexpr auto getFilledInIsNonzeroInfo(
    const IsNonzeroInfo<dim, dim> is_nonzero_original_matrix) {
  IsNonzeroInfo<dim, dim> is_nonzero;
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      is_nonzero[i][j] = is_nonzero_original_matrix[i][j];
      for (std::size_t k = 0; k < j; ++k) {
        if (is_nonzero[i][k] && is_nonzero[j][k]) {
          is_nonzero[i][j] = true;
        }
      }
    }
  }
  return is_nonzero;
}

}  // namespace ctldl
