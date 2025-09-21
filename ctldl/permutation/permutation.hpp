#pragma once

#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <std::size_t dim>
struct Permutation {
  static constexpr std::size_t size() { return dim; }

  constexpr Permutation()
      : indices([] {
          std::array<std::size_t, dim> permutation;
          fixInitIfZeroLengthArray(permutation);
          std::iota(permutation.begin(), permutation.end(), std::size_t{0});
          return permutation;
        }()) {}

  constexpr explicit Permutation(
      const std::array<std::size_t, dim>& permutation)
      : indices(permutation) {}

  constexpr Permutation(const PermutationIdentity) : Permutation() {}

  constexpr auto operator[](const std::size_t i) const { return indices[i]; }

  constexpr auto begin() const { return indices.begin(); }
  constexpr auto end() const { return indices.end(); }
  constexpr auto cbegin() const { return indices.cbegin(); }
  constexpr auto cend() const { return indices.cend(); }

  std::array<std::size_t, dim> indices;
};

template <std::size_t dim>
constexpr bool operator==(const Permutation<dim>& lhs,
                          const Permutation<dim>& rhs) {
  return std::ranges::equal(lhs, rhs);
}

}  // namespace ctldl
