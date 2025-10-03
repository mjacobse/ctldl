#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <cstddef>
#include <vector>

namespace ctldl {

constexpr auto invertPermutation(const PermutationView permutation) {
  std::vector<std::size_t> inverse_permutation(permutation.size());
  for (std::size_t i = 0; i < permutation.size(); ++i) {
    inverse_permutation[permutation[i]] = i;
  }
  return PermutationDynamic(inverse_permutation);
}

}  // namespace ctldl
