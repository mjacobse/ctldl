#pragma once

#include <ctldl/factor_data/sparsity_to_factorize_link.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_outer.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal.hpp>

#include <cstddef>

namespace ctldl {

/**
 * Represents the sparsity of a symmetric block-tridiagonal arrowhead matrix
 * with repeating sparsity and connecting link blocks between the tridiagonal
 * and outer part
 *
 * [A  B'             E']
 * [B  A  B'          E']
 * [   :  :  :        : ]
 * [      B  A  B'    E']
 * [         B  A  C' E']
 * [            C  D  F']
 * [E  E ... E  E  F  G ]
 *
 * by storing the sparsity of its repeating tridiagonal blocks A and B, its
 * linking blocks C, D and F and its outer blocks E and G.
 */
template <std::size_t dim_tridiag_, std::size_t dim_link_,
          std::size_t dim_outer_, std::size_t nnz_tridiag_diag,
          std::size_t nnz_tridiag_subdiag, std::size_t nnz_link_tridiag,
          std::size_t nnz_link_diag, std::size_t nnz_link_outer,
          std::size_t nnz_outer_subdiagonal, std::size_t nnz_outer_diagonal>
struct SparsityToFactorizeTridiagonalArrowheadLinked {
  static constexpr auto dim_tridiag = dim_tridiag_;
  static constexpr auto dim_link = dim_link_;
  static constexpr auto dim_outer = dim_outer_;
  SparsityToFactorizeTridiagonal<dim_tridiag_, nnz_tridiag_diag,
                                 nnz_tridiag_subdiag>
      tridiag;
  SparsityToFactorizeLink<dim_tridiag_, dim_link_, dim_outer_, nnz_link_tridiag,
                          nnz_link_diag, nnz_link_outer>
      link;
  SparsityToFactorizeOuter<dim_tridiag_, dim_outer_, nnz_outer_subdiagonal,
                           nnz_outer_diagonal>
      outer;
};

// needed for clang < 17 which does not do the CTAD for aggregates otherwise
template <std::size_t dim_tridiag, std::size_t dim_link, std::size_t dim_outer,
          std::size_t nnz_tridiag_diag, std::size_t nnz_tridiag_subdiag,
          std::size_t nnz_link_tridiag, std::size_t nnz_link_diag,
          std::size_t nnz_link_outer, std::size_t nnz_outer_subdiagonal,
          std::size_t nnz_outer_diagonal>
SparsityToFactorizeTridiagonalArrowheadLinked(
    SparsityToFactorizeTridiagonal<dim_tridiag, nnz_tridiag_diag,
                                   nnz_tridiag_subdiag>,
    SparsityToFactorizeLink<dim_tridiag, dim_link, dim_outer, nnz_link_tridiag,
                            nnz_link_diag, nnz_link_outer>,
    SparsityToFactorizeOuter<dim_tridiag, dim_outer, nnz_outer_subdiagonal,
                             nnz_outer_diagonal>)
    -> SparsityToFactorizeTridiagonalArrowheadLinked<
        dim_tridiag, dim_link, dim_outer, nnz_tridiag_diag, nnz_tridiag_subdiag,
        nnz_link_tridiag, nnz_link_diag, nnz_link_outer, nnz_outer_subdiagonal,
        nnz_outer_diagonal>;

}  // namespace ctldl
