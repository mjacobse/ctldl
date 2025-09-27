#pragma once

#include <ctldl/matrix/multiply.hpp>
#include <ctldl/matrix/multiply_kind.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>

namespace ctldl {

template <class MatricesA, class MatricesB, class Solution, class Rhs>
void multiplyRepeatingBlockTridiagonal(const MatricesA& matrices_A,
                                       const MatricesB& matrices_B,
                                       const Solution& solution, Rhs& rhs) {
  const auto num_repetitions = std::size_t{matrices_B.size()};
  pre(num_repetitions < std::numeric_limits<std::size_t>::max());
  pre(matrices_A.size() == num_repetitions + 1);
  pre(solution.size() == num_repetitions + 1);
  pre(rhs.size() == num_repetitions + 1);

  using enum MultiplyKind;
  multiply<Symmetric>(matrices_A[0], solution[0], rhs[0]);
  for (std::size_t i = 0; i < num_repetitions; ++i) {
    multiply<Transposed>(matrices_B[i], solution[i + 1], rhs[i]);
    multiply<Normal>(matrices_B[i], solution[i], rhs[i + 1]);
    multiply<Symmetric>(matrices_A[i + 1], solution[i + 1], rhs[i + 1]);
  }
}

}  // namespace ctldl
