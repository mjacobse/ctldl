#pragma once

#include <cstddef>

namespace ctldl {

template <class TestMatrixDistribution, class ValueGenerator>
auto generateMatrices(TestMatrixDistribution, ValueGenerator& value_generator,
                      const std::size_t num_matrices) {
  return TestMatrixDistribution::generate(value_generator, num_matrices);
}

template <class TestMatrixDistribution, class ValueGenerator>
auto generateMatrix(TestMatrixDistribution distribution,
                    ValueGenerator& value_generator) {
  return generateMatrices(distribution, value_generator, 1)[0];
}

}  // namespace ctldl
