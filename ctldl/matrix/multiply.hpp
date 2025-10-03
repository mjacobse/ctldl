#pragma once

#include <ctldl/matrix/multiply_kind.hpp>

#include <cstddef>
#include <type_traits>

namespace ctldl {

template <MultiplyKind kind, class Matrix, class Solution, class Rhs>
void multiply(const Matrix& matrix, const Solution& solution, Rhs&& rhs) {
  using ValueRhs = std::remove_reference_t<decltype(rhs[0])>;

  constexpr auto nnz = std::size_t{Matrix::sparsity.entries().size()};
  for (std::size_t entry_index = 0; entry_index < nnz; ++entry_index) {
    const auto row_index = Matrix::sparsity.entries()[entry_index].row_index;
    const auto col_index = Matrix::sparsity.entries()[entry_index].col_index;
    const auto value = matrix.valueAt(entry_index);
    if constexpr (kind & MultiplyKind::Normal) {
      rhs[row_index] += static_cast<ValueRhs>(value) *
                        static_cast<ValueRhs>(solution[col_index]);
    }
    if constexpr (kind == MultiplyKind::Symmetric) {
      if (row_index == col_index) {
        continue;
      }
    }
    if constexpr (kind & MultiplyKind::Transposed) {
      rhs[col_index] += static_cast<ValueRhs>(value) *
                        static_cast<ValueRhs>(solution[row_index]);
    }
  }
}

}  // namespace ctldl
