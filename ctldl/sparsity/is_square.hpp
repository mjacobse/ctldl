#pragma once

namespace ctldl {

template <class Sparsity>
constexpr bool isSquare() {
  return Sparsity::numRows() == Sparsity::numCols();
}

template <class Sparsity>
constexpr bool isSquare(const Sparsity& /*sparsity*/) {
  return isSquare<Sparsity>();
}

}  // namespace ctldl
