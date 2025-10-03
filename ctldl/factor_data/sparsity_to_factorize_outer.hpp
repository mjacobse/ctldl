#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/ctad.hpp>

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
template <std::size_t dim_inner_, std::size_t dim_outer,
          std::size_t nnz_subdiag, std::size_t nnz_diag>
struct SparsityToFactorizeOuter {
  static constexpr auto dim = dim_outer;
  static constexpr auto dim_inner = dim_inner_;
  Sparsity<nnz_subdiag, dim_outer, dim_inner_> subdiag;
  Sparsity<nnz_diag, dim_outer, dim_outer> diag;
  Permutation<dim_outer> permutation = PermutationIdentity{};
};

template <class Subdiag, class Diag, class PermutationIn = PermutationIdentity>
SparsityToFactorizeOuter(Subdiag, Diag, PermutationIn = PermutationIdentity{})
    -> SparsityToFactorizeOuter<
        ctad_t<Sparsity, Subdiag>::numCols(), ctad_t<Sparsity, Diag>::numRows(),
        ctad_t<Sparsity, Subdiag>::nnz(), ctad_t<Sparsity, Diag>::nnz()>;

template <std::size_t dim_inner, std::size_t dim_outer>
constexpr auto makeEmptySparsityToFactorizeOuter() {
  return SparsityToFactorizeOuter{makeEmptySparsity<dim_outer, dim_inner>(),
                                  makeEmptySparsity<dim_outer, dim_outer>()};
}

}  // namespace ctldl
