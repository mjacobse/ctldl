#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/ctad.hpp>

#include <cstddef>

namespace ctldl {

/**
 * Represents the sparsity of a symmetric block-tridiagonal matrix with
 * repeating sparsity
 *
 * [A  B'         ]
 * [B  A  B'      ]
 * [   :  :  :    ]
 * [      B  A  B']
 * [         B  A ]
 *
 * by storing the sparsity of its repeating diagonal blocks A and subdiagonal
 * blocks B.
 *
 * Additionally, a freely choosable symmetric permutation (identical for each
 * block) for efficient factorization is included as well.
 */
template <std::size_t dim_, std::size_t nnz_diagonal,
          std::size_t nnz_subdiagonal>
struct SparsityToFactorizeTridiagonal {
  static constexpr auto dim = dim_;
  Sparsity<nnz_diagonal, dim, dim> diag;
  Sparsity<nnz_subdiagonal, dim, dim> subdiag;
  PermutationStatic<dim> permutation = PermutationIdentity{};
};

template <class SparsityDiag, class SparsitySubdiag,
          class PermutationIn = PermutationIdentity>
SparsityToFactorizeTridiagonal(SparsityDiag, SparsitySubdiag,
                               PermutationIn = PermutationIdentity{})
    -> SparsityToFactorizeTridiagonal<ctad_t<Sparsity, SparsityDiag>::numRows(),
                                      ctad_t<Sparsity, SparsityDiag>::nnz(),
                                      ctad_t<Sparsity, SparsitySubdiag>::nnz()>;

}  // namespace ctldl
