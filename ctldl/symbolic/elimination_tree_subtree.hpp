#pragma once

#include <ctldl/symbolic/elimination_tree.hpp>
#include <ctldl/utility/contracts.hpp>

#include <algorithm>
#include <cstddef>


namespace ctldl {

constexpr EliminationTree subtreeBack(
    const EliminationTree& tree, const std::size_t dim_subtree) {
  pre(dim_subtree <= tree.parent.size());
  const auto dim = std::size_t{tree.parent.size()};
  EliminationTree subtree(dim_subtree);
  std::transform(tree.parent.cend() - static_cast<std::ptrdiff_t>(dim_subtree),
                 tree.parent.cend(), subtree.parent.begin(),
                 [&](const std::size_t node) {
                   return node == tree.no_parent ? subtree.no_parent
                                                 : node - (dim - dim_subtree);
                 });
  return subtree;
}

constexpr void clearSubtreeBack(EliminationTree& tree,
                                const std::size_t dim_subtree) {
  std::fill(tree.parent.end() - static_cast<std::ptrdiff_t>(dim_subtree),
            tree.parent.end(), tree.no_parent);
}

constexpr void setSubtreeFront(EliminationTree& tree,
                               const EliminationTree& subtree) {
  pre(subtree.parent.size() <= tree.parent.size());
  std::ranges::copy(subtree.parent, tree.parent.begin());
}

}  // namespace ctldl
