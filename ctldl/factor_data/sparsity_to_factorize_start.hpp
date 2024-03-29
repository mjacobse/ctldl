#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/ctad.hpp>

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
template <std::size_t dim_start, std::size_t dim_next_, std::size_t dim_outer_,
          std::size_t nnz_start, std::size_t nnz_next, std::size_t nnz_outer>
struct SparsityToFactorizeStart {
  static constexpr auto dim = dim_start;
  static constexpr auto dim_next = dim_next_;
  static constexpr auto dim_outer = dim_outer_;
  Sparsity<nnz_start, dim_start, dim_start> diag;
  Sparsity<nnz_next, dim_next_, dim_start> next;
  Sparsity<nnz_outer, dim_outer_, dim_start> outer;
  Permutation<dim_start> permutation = PermutationIdentity{};
};

template <class Diag, class Next, class Outer,
          class PermutationIn = PermutationIdentity>
SparsityToFactorizeStart(Diag, Next, Outer,
                         PermutationIn = PermutationIdentity{})
    -> SparsityToFactorizeStart<
        ctad_t<Sparsity, Diag>::num_rows, ctad_t<Sparsity, Next>::num_rows,
        ctad_t<Sparsity, Outer>::num_rows, ctad_t<Sparsity, Diag>::nnz,
        ctad_t<Sparsity, Next>::nnz, ctad_t<Sparsity, Outer>::nnz>;

template <std::size_t dim_start, std::size_t dim_next, std::size_t dim_outer>
constexpr auto makeEmptySparsityToFactorizeStart() {
  return SparsityToFactorizeStart{makeEmptySparsity<dim_start, dim_start>(),
                                  makeEmptySparsity<dim_next, dim_start>(),
                                  makeEmptySparsity<dim_outer, dim_start>()};
}

}  // namespace ctldl
