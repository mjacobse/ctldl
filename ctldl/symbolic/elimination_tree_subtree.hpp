#pragma once

#include <ctldl/symbolic/elimination_tree.hpp>

#include <algorithm>
#include <cstddef>


namespace ctldl {

template <std::size_t dim_subtree, std::size_t dim>
constexpr EliminationTree<dim_subtree> subtreeBack(
    const EliminationTree<dim>& tree) {
  static_assert(dim_subtree <= dim);
  EliminationTree<dim_subtree> subtree;
  std::transform(tree.parent.cend() - dim_subtree, tree.parent.cend(),
                 subtree.parent.begin(), [&](const std::size_t node) {
                   return node == tree.no_parent ? subtree.no_parent
                                                 : node - (dim - dim_subtree);
                 });
  return subtree;
}

template <std::size_t dim>
constexpr void clearSubtreeBack(EliminationTree<dim>& tree,
                                const std::size_t dim_subtree) {
  std::fill(tree.parent.end() - dim_subtree, tree.parent.end(), tree.no_parent);
}

template <std::size_t dim, std::size_t dim_subtree>
constexpr void setSubtreeFront(EliminationTree<dim>& tree,
                               const EliminationTree<dim_subtree>& subtree) {
  static_assert(dim_subtree <= dim);
  std::copy(subtree.parent.cbegin(), subtree.parent.cend(),
            tree.parent.begin());
}

}  // namespace ctldl
