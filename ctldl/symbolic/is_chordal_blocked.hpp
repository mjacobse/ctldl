#pragma once

#include <ctldl/sparsity/is_lower_triangle.hpp>
#include <ctldl/sparsity/is_sparsity_equal.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/filled_in_sparsity_blocked.hpp>

namespace ctldl {

template <class Sparsity11, class Sparsity21, class Sparsity22>
constexpr bool isChordalBlocked(const Sparsity11& sparsity11,
                                const Sparsity21& sparsity21,
                                const Sparsity22& sparsity22) {
  static_assert(isSquare<Sparsity11>());
  static_assert(isSquare<Sparsity22>());
  static_assert(Sparsity21::num_cols == Sparsity11::num_cols);
  static_assert(Sparsity22::num_rows == Sparsity21::num_rows);
  if (!isLowerTriangle(sparsity11) || !isLowerTriangle(sparsity22)) {
    return false;
  }
  const auto nnz_filled =
      getFilledInNumNonZerosBlocked(sparsity11, sparsity21, sparsity22);
  return nnz_filled.block11 == sparsity11.nnz &&
         nnz_filled.block21 == sparsity21.nnz &&
         nnz_filled.block22 == sparsity22.nnz;
};

}  // namespace ctldl
