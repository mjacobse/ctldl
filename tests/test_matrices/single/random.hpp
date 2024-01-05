#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <vector>

namespace ctldl {

template <auto sparsity_in, class Value>
struct TestMatrixRandom {
  static constexpr const char* description() {
    return "Random";
  }

  struct Matrix {
    static constexpr auto sparsity = Sparsity(sparsity_in);
    std::array<Value, sparsity.nnz> values;
    Value valueAt(const std::size_t entry_index) const {
      return values[entry_index];
    }
  };

  template <class ValueGenerator>
  static auto generate(ValueGenerator& value_generator,
                       const std::size_t num_matrices) {
    std::vector<Matrix> matrices(num_matrices);
    const auto standard_deviation = Value{1.0};
    std::normal_distribution<> distribution(0.0, standard_deviation);
    for (auto& matrix : matrices) {
      std::generate(matrix.values.begin(), matrix.values.end(),
                    [&distribution, &value_generator] {
                      return distribution(value_generator);
                    });
    }
    return matrices;
  }
};

}  // namespace ctldl
