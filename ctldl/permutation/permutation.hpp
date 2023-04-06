#pragma once

#include <ctldl/permutation/permutation_identity.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>

namespace ctldl {

template <std::size_t dim>
class Permutation {
 public:
  constexpr Permutation()
      : m_permutation([] {
          std::array<std::size_t, dim> permutation;
          std::iota(permutation.begin(), permutation.end(), std::size_t{0});
          return permutation;
        }()) {}

  constexpr explicit Permutation(
      const std::array<std::size_t, dim>& permutation)
      : m_permutation(permutation) {}

  constexpr Permutation(const PermutationIdentity::ContructTag)
      : Permutation() {}

  constexpr auto operator[](const std::size_t i) const {
    return m_permutation[i];
  }

  constexpr auto begin() const { return m_permutation.begin(); }
  constexpr auto end() const { return m_permutation.end(); }

 private:
  std::array<std::size_t, dim> m_permutation;
};

template <std::size_t dim>
constexpr bool operator==(const Permutation<dim>& lhs,
                          const Permutation<dim>& rhs) {
  return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

}  // namespace ctldl
