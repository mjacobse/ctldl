#pragma once

#include <ctldl/matrix/empty_matrix_input.hpp>

#include <cstddef>

namespace ctldl {

template <class Diag, class Next, class Outer>
struct MatrixStart {
  Diag diag;
  Next next;
  Outer outer;
};

template <std::size_t dim_start, std::size_t dim_next, std::size_t dim_outer>
constexpr auto makeEmptyMatrixStart() {
  return MatrixStart{EmptyMatrixInput<dim_start, dim_start>{},
                     EmptyMatrixInput<dim_next, dim_start>{},
                     EmptyMatrixInput<dim_outer, dim_start>{}};
}

}  // namespace ctldl
