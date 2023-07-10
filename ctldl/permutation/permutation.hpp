#pragma once

#include <ctldl/permutation/permutation_identity.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <std::size_t dim>
struct Permutation {
  static constexpr std::size_t size() { return dim; }

  constexpr explicit Permutation()
      : indices([] {
          std::array<std::size_t, dim> permutation;
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

  std::array<std::size_t, dim> indices;
};

template <std::size_t dim>
constexpr bool operator==(const Permutation<dim>& lhs,
                          const Permutation<dim>& rhs) {
  return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

}  // namespace ctldl
