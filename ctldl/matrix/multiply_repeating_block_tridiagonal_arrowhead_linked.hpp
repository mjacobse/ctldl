#pragma once

#include <ctldl/matrix/multiply.hpp>
#include <ctldl/matrix/multiply_kind.hpp>
#include <ctldl/matrix/multiply_repeating_block_tridiagonal.hpp>

#include <cstddef>

namespace ctldl {

template <class Matrix, class Solution, class Rhs>
void multiplyRepeatingBlockTridiagonalArrowheadLinked(const Matrix& matrix,
                                                      const Solution& solution,
                                                      Rhs& rhs) {
  using enum MultiplyKind;

  multiply<Symmetric>(matrix.start.diag, solution.start, rhs.start);
  multiply<Transposed>(matrix.start.next, solution.tridiag.front(), rhs.start);
  multiply<Transposed>(matrix.start.outer, solution.outer, rhs.start);

  multiply<Normal>(matrix.start.next, solution.start, rhs.tridiag.front());
  multiplyRepeatingBlockTridiagonal(matrix.tridiag.diag, matrix.tridiag.subdiag,
                                    solution.tridiag, rhs.tridiag);

  const auto num_repetitions = std::size_t{matrix.tridiag.subdiag.size()};
  for (std::size_t i = 0; i <= num_repetitions; ++i) {
    multiply<Transposed>(matrix.outer.subdiag[i], solution.outer,
                         rhs.tridiag[i]);
    multiply<Normal>(matrix.outer.subdiag[i], solution.tridiag[i], rhs.outer);
  }
  multiply<Transposed>(matrix.link.prev, solution.link, rhs.tridiag.back());

  multiply<Normal>(matrix.link.prev, solution.tridiag.back(), rhs.link);
  multiply<Transposed>(matrix.link.next, solution.outer, rhs.link);
  multiply<Symmetric>(matrix.link.diag, solution.link, rhs.link);

  multiply<Normal>(matrix.start.outer, solution.start, rhs.outer);
  multiply<Normal>(matrix.link.next, solution.link, rhs.outer);
  multiply<Symmetric>(matrix.outer.diag, solution.outer, rhs.outer);
}

}  // namespace ctldl
