#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/sparsity.hpp>
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
template <std::size_t dim_prev_, std::size_t dim_link, std::size_t dim_next_,
          std::size_t nnz_prev, std::size_t nnz_link, std::size_t nnz_next>
struct SparsityToFactorizeLink {
  static constexpr auto dim_prev = dim_prev_;
  static constexpr auto dim = dim_link;
  static constexpr auto dim_next = dim_next_;
  SparsityStatic<nnz_prev, dim_link, dim_prev_> prev;
  SparsityStatic<nnz_link, dim_link, dim_link> diag;
  SparsityStatic<nnz_next, dim_next_, dim_link> next;
  PermutationStatic<dim_link> permutation = PermutationIdentity{};
};

template <class Prev, class Diag, class Next,
          class PermutationIn = PermutationIdentity>
SparsityToFactorizeLink(Prev, Diag, Next, PermutationIn = PermutationIdentity{})
    -> SparsityToFactorizeLink<ctad_t<SparsityStatic, Prev>::numCols(),
                               ctad_t<SparsityStatic, Diag>::numRows(),
                               ctad_t<SparsityStatic, Next>::numRows(),
                               ctad_t<SparsityStatic, Prev>::nnz(),
                               ctad_t<SparsityStatic, Diag>::nnz(),
                               ctad_t<SparsityStatic, Next>::nnz()>;

template <std::size_t dim_prev, std::size_t dim_link, std::size_t dim_next>
constexpr auto makeEmptySparsityToFactorizeLink() {
  return SparsityToFactorizeLink{makeEmptySparsityStatic<dim_link, dim_prev>(),
                                 makeEmptySparsityStatic<dim_link, dim_link>(),
                                 makeEmptySparsityStatic<dim_next, dim_link>()};
}

}  // namespace ctldl
