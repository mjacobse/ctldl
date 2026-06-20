#pragma once

#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <meta>
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

  constexpr explicit PermutationView(const std::span<const std::size_t> indices)
      : m_indices(indices) {
    // TODO: check invariant
  }

  friend class PermutationViewStructural;
};

constexpr bool operator==(const PermutationView lhs,
                          const PermutationView rhs) {
  return std::ranges::equal(lhs, rhs);
}

class PermutationViewStructural {
 public:
  const std::size_t* m_indices_do_not_touch;
  std::size_t m_size_do_not_touch;

  consteval explicit(false)
      PermutationViewStructural(const PermutationView permutation)
      : m_indices_do_not_touch(std::define_static_array(permutation).data()),
        m_size_do_not_touch(permutation.size()) {}

  consteval explicit(false)
      PermutationViewStructural(const PermutationDynamic& permutation)
      : PermutationViewStructural(PermutationView(permutation)) {}

  template <std::size_t dim>
  consteval explicit(false)
      PermutationViewStructural(const PermutationStatic<dim>& permutation)
      : PermutationViewStructural(PermutationView(permutation)) {}

  constexpr std::size_t size() const { return m_size_do_not_touch; }

  constexpr auto operator[](const std::size_t i) const {
    return m_indices_do_not_touch[i];
  }

  constexpr auto begin() const { return m_indices_do_not_touch; }
  constexpr auto end() const { return m_indices_do_not_touch + size(); }

  constexpr operator PermutationView() const {
    return PermutationView(std::span<const std::size_t>(begin(), end()));
  }
};

}  // namespace ctldl
