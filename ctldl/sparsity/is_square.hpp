#pragma once

#include <ctldl/sparsity/sparsity.hpp>

namespace ctldl {

constexpr bool isSquare(const SparsityView sparsity) {
  return sparsity.numRows() == sparsity.numCols();
}

}  // namespace ctldl
