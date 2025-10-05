#pragma once

#include <ctldl/sparsity/is_lower_triangle.hpp>
#include <ctldl/sparsity/is_sparsity_equal.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/filled_in_sparsity_blocked.hpp>
#include <ctldl/utility/contracts.hpp>

namespace ctldl {

template <class Sparsity11, class Sparsity21, class Sparsity22>
constexpr bool isChordalBlocked(const Sparsity11& sparsity11,
                                const Sparsity21& sparsity21,
                                const Sparsity22& sparsity22) {
  pre(isSquare(sparsity11));
  pre(isSquare(sparsity22));
  static_assert(Sparsity21::numCols() == Sparsity11::numCols());
  static_assert(Sparsity22::numRows() == Sparsity21::numRows());
  if (!isLowerTriangle(sparsity11) || !isLowerTriangle(sparsity22)) {
    return false;
  }
  const auto nnz_filled =
      getFilledInNumNonZerosBlocked(sparsity11, sparsity21, sparsity22);
  return nnz_filled.block11 == sparsity11.nnz() &&
         nnz_filled.block21 == sparsity21.nnz() &&
         nnz_filled.block22 == sparsity22.nnz();
};

}  // namespace ctldl
