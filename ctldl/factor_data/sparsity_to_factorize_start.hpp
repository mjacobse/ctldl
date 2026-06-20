#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/all_equal.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>

namespace ctldl {

/**
 * Represents the start of a next and an outer part
 *
 * [A  B'             C']
 * [B  *  *' *' *' *' *']
 * [   *  *  *' *' *' *']
 * [   *  *  *  *' *' *']
 * [   *  *  *  *  *' *']
 * [   *  *  *  *  *  *']
 * [C  *  *  *  *  *  *']
 *
 * by storing the sparsity of its diagonal block A, its connecting block to the
 * the next part B and its connecting block to the outer part C.
 *
 * Additionally, a freely choosable permutation (symmetric permutation of A,
 * row permutation of B, row permutation of C) for efficient factorization is
 * included as well.
 */
struct SparsityToFactorizeStart {
  SparsityViewStructural diag;
  SparsityViewStructural next;
  SparsityViewStructural outer;
  PermutationViewStructural permutation;

  constexpr bool has_consistent_dim() const {
    return all_equal({diag.numRows(), diag.numCols(), next.numCols(),
                      outer.numCols(), permutation.size()});
  }

  constexpr std::size_t dim() const {
    pre(has_consistent_dim());
    return diag.numCols();
  }

  constexpr std::size_t dim_next() const { return next.numRows(); }
  constexpr std::size_t dim_outer() const { return outer.numRows(); }
};

template <std::size_t dim_start, std::size_t dim_next, std::size_t dim_outer>
constexpr auto makeEmptySparsityToFactorizeStart() {
  return SparsityToFactorizeStart{
      makeEmptySparsityStatic<dim_start, dim_start>(),
      makeEmptySparsityStatic<dim_next, dim_start>(),
      makeEmptySparsityStatic<dim_outer, dim_start>(),
      PermutationDynamic(dim_start)};
}

}  // namespace ctldl
