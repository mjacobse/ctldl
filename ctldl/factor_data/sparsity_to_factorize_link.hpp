#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/all_equal.hpp>
#include <ctldl/utility/contracts.hpp>
#include <ctldl/utility/ctad.hpp>

#include <cstddef>

namespace ctldl {

/**
 * Represents the link between a previous and a next part
 *
 * [*  *' *'    *' *' *']
 * [*  *  *'    *' *' *']
 * [*  *  *  A' *' *' *']
 * [      A  B  C'      ]
 * [*  *  *  C  *  *' *']
 * [*  *  *     *  *  *']
 * [*  *  *     *  *  *']
 *
 * by storing the sparsity of its connecting block to the previous part A,
 * its own diagonal block B and its connecting block to the outer part C.
 *
 * Additionally, a freely choosable permutation (symmetric permutation of B,
 * row permutation of A, column permutation of C) for efficient factorization is
 * included as well.
 */
struct SparsityToFactorizeLink {
  SparsityViewStructural prev;
  SparsityViewStructural diag;
  SparsityViewStructural next;
  PermutationViewStructural permutation;

  constexpr bool has_consistent_dim() const {
    return all_equal({diag.numRows(), diag.numCols(), prev.numRows(),
                      next.numCols(), permutation.size()});
  }

  constexpr std::size_t dim() const {
    pre(has_consistent_dim());
    return diag.numRows();
  }

  constexpr std::size_t dim_prev() const { return prev.numCols(); }
  constexpr std::size_t dim_next() const { return next.numRows(); }
};

template <std::size_t dim_prev, std::size_t dim_link, std::size_t dim_next>
constexpr auto makeEmptySparsityToFactorizeLink() {
  return SparsityToFactorizeLink{makeEmptySparsityStatic<dim_link, dim_prev>(),
                                 makeEmptySparsityStatic<dim_link, dim_link>(),
                                 makeEmptySparsityStatic<dim_next, dim_link>(),
                                 PermutationDynamic(dim_link)};
}

}  // namespace ctldl
