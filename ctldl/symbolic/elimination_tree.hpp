#pragma once

#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>


namespace ctldl {

template <std::size_t dim>
struct EliminationTree {
  static constexpr auto no_parent = std::numeric_limits<std::size_t>::max();

  constexpr EliminationTree() {
    fixInitIfZeroLengthArray(parent);
    std::ranges::fill(parent, no_parent);
  };

  constexpr bool hasParent(const std::size_t node) const {
    return parent[node] != no_parent;
  }

  constexpr auto getParent(const std::size_t node) const {
    assert(hasParent(node));
    return parent[node];
  }

  std::array<std::size_t, dim> parent;
};

template <std::size_t dim>
constexpr bool operator==(const EliminationTree<dim>& lhs,
                          const EliminationTree<dim>& rhs) {
  return std::ranges::equal(lhs.parent, rhs.parent);
}

}  // namespace ctldl
