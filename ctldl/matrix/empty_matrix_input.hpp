#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>

#include <cassert>
#include <cstddef>

namespace ctldl {

template <std::size_t num_rows, std::size_t num_cols>
struct EmptyMatrixInput {
  static constexpr auto sparsity =
      makeEmptySparsityStaticCSR<num_rows, num_cols>();

  constexpr EmptyMatrixInput() = default;

  template <class Value>
  constexpr explicit EmptyMatrixInput(const std::array<Value, 0> /*values*/) {}

  constexpr double valueAt(const std::size_t /*i*/) const {
    assert(false);
    return 0.0;
  }
};

}  // namespace ctldl
