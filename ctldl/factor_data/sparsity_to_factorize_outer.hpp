#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/all_equal.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>

namespace ctldl {

/**
 * Represents the sparsity of the outer part of a matrix
 *
 * [*  A']
 * [A  B ]
 *
 * by storing the sparsity of its subdiagonal block B and its diagonal block A.
 *
 * Additionally, a freely choosable permutation (symmetric permutation of B,
 * row permutation of A) for efficient factorization is included as well.
 *
 * This is essentially the same as SparsityToFactorizeTail, which is why we
 * just type-alias that. The name is more descriptive though, for example if
 * describing the outer part of a tridiagonal arrowhead matrix
 *
 * [*  *'          A']
 * [*  *  *'       A']
 * [   :  :  :     : ]
 * [      *  *  *' A']
 * [         *  *  A']
 * [A  A ... A  A  B ]
 *
 * with repeating A.
 */
struct SparsityToFactorizeOuter {
  SparsityViewStructural subdiag;
  SparsityViewStructural diag;
  PermutationViewStructural permutation;

  constexpr bool has_consistent_dim() const {
    return all_equal({diag.numRows(), diag.numCols(), subdiag.numRows(),
                      permutation.size()});
  }

  constexpr std::size_t dim() const {
    pre(has_consistent_dim());
    return diag.numRows();
  }

  constexpr std::size_t dim_inner() const { return subdiag.numCols(); }
};

template <std::size_t dim_inner, std::size_t dim_outer>
constexpr auto makeEmptySparsityToFactorizeOuter() {
  return SparsityToFactorizeOuter{
      makeEmptySparsityStatic<dim_outer, dim_inner>(),
      makeEmptySparsityStatic<dim_outer, dim_outer>(),
      PermutationDynamic(dim_outer)};
}

}  // namespace ctldl
