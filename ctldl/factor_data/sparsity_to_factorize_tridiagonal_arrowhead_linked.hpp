#pragma once

#include <ctldl/factor_data/sparsity_to_factorize_link.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_outer.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_start.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal.hpp>
#include <ctldl/utility/all_equal.hpp>
#include <ctldl/utility/contracts.hpp>

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
struct SparsityToFactorizeTridiagonalArrowheadLinked {
  SparsityToFactorizeStart start;
  SparsityToFactorizeTridiagonal tridiag;
  SparsityToFactorizeLink link;
  SparsityToFactorizeOuter outer;

  constexpr bool has_consistent_dim_tridiag() const {
    return tridiag.has_consistent_dim() &&
           all_equal({start.dim_next(), tridiag.dim(), link.dim_prev(),
                      outer.dim_inner()});
  }

  constexpr bool has_consistent_dim_outer() const {
    return outer.has_consistent_dim() &&
           all_equal({start.dim_outer(), link.dim_next(), outer.dim()});
  }

  constexpr std::size_t dim_start() const { return start.dim(); }

  constexpr std::size_t dim_tridiag() const {
    pre(has_consistent_dim_tridiag());
    return tridiag.dim();
  }

  constexpr std::size_t dim_link() const { return link.dim(); }

  constexpr std::size_t dim_outer() const {
    pre(has_consistent_dim_outer());
    return outer.dim();
  }
};

}  // namespace ctldl
