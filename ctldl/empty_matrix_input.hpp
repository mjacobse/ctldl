#pragma once

#include <ctldl/sparsity/sparsity.hpp>

#include <cstddef>

namespace ctldl {

template <std::size_t num_rows, std::size_t num_cols>
struct EmptyMatrixInput {
  static constexpr auto sparsity = makeEmptySparsity<num_rows, num_cols>();
};

}  // namespace ctldl
