#pragma once

#include <ctldl/sparsity/is_lower_triangle.hpp>
#include <ctldl/sparsity/is_sparsity_equal.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/filled_in_sparsity_blocked.hpp>
#include <ctldl/utility/contracts.hpp>

namespace ctldl {

constexpr bool isChordalBlocked(const SparsityViewCSR sparsity11,
                                const SparsityViewCSR sparsity21,
                                const SparsityViewCSR sparsity22) {
  pre(isSquare(sparsity11));
  pre(isSquare(sparsity22));
  pre(sparsity21.numCols() == sparsity11.numCols());
  pre(sparsity22.numRows() == sparsity21.numRows());
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
