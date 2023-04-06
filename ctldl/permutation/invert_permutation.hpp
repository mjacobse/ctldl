#pragma once

#include <ctldl/permutation/permutation.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <std::size_t dim>
constexpr auto invertPermutation(const Permutation<dim>& permutation) {
  std::array<std::size_t, dim> inverse_permutation;
  for (std::size_t i = 0; i < dim; ++i) {
    inverse_permutation[permutation[i]] = i;
  }
  return Permutation<dim>{inverse_permutation};
}

}  // namespace ctldl
