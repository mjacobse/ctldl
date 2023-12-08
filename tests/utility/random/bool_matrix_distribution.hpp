#pragma once

#include "tests/utility/random/bernoulli_distribution.hpp"

#include <ctldl/sparsity/bool_matrix.hpp>

#include <cstddef>

namespace ctldl {

template <std::size_t num_rows, std::size_t num_cols>
struct BoolMatrixDistribution {
  using result_type = BoolMatrix<num_rows, num_cols>;

  template <class Generator>
  constexpr auto operator()(Generator& generator) const {
    BoolMatrix<num_rows, num_cols> is_nonzero;
    BernoulliDistribution distribution(0.5);
    for (std::size_t i = 0; i < num_rows; ++i) {
      for (std::size_t j = 0; j < num_cols; ++j) {
        is_nonzero.values[i][j] = distribution(generator);
      }
    }
    return is_nonzero;
  }
};

}  // namespace ctldl
