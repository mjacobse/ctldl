#pragma once

#include <ctldl/matrix/empty_matrix_input.hpp>
#include <ctldl/utility/dummy_range.hpp>

#include <cstddef>

namespace ctldl {

template <class Subdiag, class Diag>
struct MatrixOuter {
  Subdiag subdiag;
  Diag diag;
};

template <std::size_t num_rows, std::size_t num_cols>
constexpr auto makeEmptyMatrixRepeatingOuter(
    const std::size_t num_repetitions) {
  DummyRange<EmptyMatrixInput<num_rows, num_cols>> dummy_subdiag(
      num_repetitions + 1);
  EmptyMatrixInput<num_rows, num_rows> dummy_diag;
  return MatrixOuter{dummy_subdiag, dummy_diag};
}

}  // namespace ctldl
