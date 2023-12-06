#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class PermutationIn>
constexpr auto invertPermutation(const PermutationIn& permutation) {
  using std::size;
  constexpr auto dim = std::size_t{size(PermutationIn{})};

  std::array<std::size_t, dim> inverse_permutation;
  fixInitIfZeroLengthArray(inverse_permutation);
  for (std::size_t i = 0; i < dim; ++i) {
    inverse_permutation[permutation[i]] = i;
  }
  return Permutation<dim>{inverse_permutation};
}

}  // namespace ctldl
