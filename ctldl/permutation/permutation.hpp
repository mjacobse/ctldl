#pragma once

#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>
#include <span>
#include <vector>

namespace ctldl {

template <std::size_t dim>
struct PermutationStatic {
  static constexpr std::size_t size() { return dim; }

  constexpr PermutationStatic()
      : m_indices_do_not_touch([] {
          std::array<std::size_t, dim> permutation;
          fixInitIfZeroLengthArray(permutation);
          std::iota(permutation.begin(), permutation.end(), std::size_t{0});
          return permutation;
        }()) {}

  constexpr explicit PermutationStatic(
      const std::array<std::size_t, dim>& permutation)
      : m_indices_do_not_touch(permutation) {}

  constexpr PermutationStatic(const PermutationIdentity)
      : PermutationStatic() {}

  constexpr auto operator[](const std::size_t i) const { return indices()[i]; }

  constexpr auto begin() const { return indices().begin(); }
  constexpr auto end() const { return indices().end(); }

  // only public to allow usage as NTTP
  std::array<std::size_t, dim> m_indices_do_not_touch;
  constexpr const std::array<std::size_t, dim>& indices() const {
    return m_indices_do_not_touch;
  }
};

class PermutationDynamic {
 public:
  constexpr explicit PermutationDynamic(
      const std::span<const std::size_t> indices)
      : m_indices(indices.begin(), indices.end()) {}

  constexpr explicit PermutationDynamic(const std::size_t dim)
      : m_indices([dim] {
          std::vector<std::size_t> permutation(dim);
          std::iota(permutation.begin(), permutation.end(), std::size_t{0});
          return permutation;
        }()) {}

  constexpr const std::vector<std::size_t>& indices() const {
    return m_indices;
  }

  constexpr std::size_t size() const { return m_indices.size(); }

  constexpr auto operator[](const std::size_t i) const { return m_indices[i]; }

  constexpr auto begin() const { return indices().begin(); }
  constexpr auto end() const { return indices().end(); }

 private:
  std::vector<std::size_t> m_indices;
};

class PermutationView {
 public:
  template <std::size_t dim>
  constexpr PermutationView(const PermutationStatic<dim>& permutation)
      : m_indices(permutation.indices()) {}

  constexpr PermutationView(const PermutationDynamic& permutation)
      : m_indices(permutation.indices()) {}

  constexpr std::size_t size() const { return m_indices.size(); }

  constexpr auto operator[](const std::size_t i) const { return m_indices[i]; }

  constexpr auto begin() const { return m_indices.begin(); }
  constexpr auto end() const { return m_indices.end(); }

 private:
  std::span<const std::size_t> m_indices;
};

constexpr bool operator==(const PermutationView lhs,
                          const PermutationView rhs) {
  return std::ranges::equal(lhs, rhs);
}

}  // namespace ctldl
