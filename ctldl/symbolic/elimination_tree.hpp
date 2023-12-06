#pragma once

#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>


namespace ctldl {

template <std::size_t dim>
struct EliminationTree {
  static constexpr auto no_parent = std::numeric_limits<std::size_t>::max();

  constexpr EliminationTree() {
    fixInitIfZeroLengthArray(parent);
    std::fill(parent.begin(), parent.end(), no_parent);
  };

  constexpr bool hasParent(const std::size_t node) const {
    return parent[node] != no_parent;
  }

  std::array<std::size_t, dim> parent;
};

template <std::size_t dim>
constexpr bool operator==(const EliminationTree<dim>& lhs,
                          const EliminationTree<dim>& rhs) {
  return std::equal(lhs.parent.cbegin(), lhs.parent.cend(),
                    rhs.parent.cbegin());
}

}  // namespace ctldl
