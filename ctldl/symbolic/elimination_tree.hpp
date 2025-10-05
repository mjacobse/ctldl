#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>


namespace ctldl {

struct EliminationTree {
  static constexpr auto no_parent = std::numeric_limits<std::size_t>::max();

  constexpr EliminationTree(const std::size_t dim) : parent(dim, no_parent) {};

  constexpr bool hasParent(const std::size_t node) const {
    return parent[node] != no_parent;
  }

  constexpr auto getParent(const std::size_t node) const {
    assert(hasParent(node));
    return parent[node];
  }

  std::vector<std::size_t> parent;
};

constexpr bool operator==(const EliminationTree& lhs,
                          const EliminationTree& rhs) {
  return std::ranges::equal(lhs.parent, rhs.parent);
}

}  // namespace ctldl
