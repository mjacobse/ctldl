#pragma once

#include <ctldl/matrix/empty_matrix_input.hpp>

#include <cstddef>

namespace ctldl {

template <class Prev, class Diag, class Next>
struct MatrixLink {
  Prev prev;
  Diag diag;
  Next next;
};

template <std::size_t dim_prev, std::size_t dim_link, std::size_t dim_next>
constexpr auto makeEmptyMatrixLink() {
  return MatrixLink{EmptyMatrixInput<dim_link, dim_prev>{},
                    EmptyMatrixInput<dim_link, dim_link>{},
                    EmptyMatrixInput<dim_next, dim_link>{}};
}

}  // namespace ctldl
