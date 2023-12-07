#pragma once

#include <cstddef>

namespace ctldl {

template <class TestMatrixDistribution>
auto generateMatrices(TestMatrixDistribution, const std::size_t num_matrices) {
  return TestMatrixDistribution::generate(num_matrices);
}

template <class TestMatrixDistribution>
auto generateMatrix(TestMatrixDistribution distribution) {
  return generateMatrices(distribution, 1)[0];
}

}  // namespace ctldl
