#pragma once

namespace ctldl {

template <class Sparsity>
constexpr bool isSquare() {
  return Sparsity::num_rows == Sparsity::num_cols;
}

template <class Sparsity>
constexpr bool isSquare(const Sparsity& /*sparsity*/) {
  return isSquare<Sparsity>();
}

}  // namespace ctldl
