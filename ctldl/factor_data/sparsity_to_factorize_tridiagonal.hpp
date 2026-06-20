#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/utility/all_equal.hpp>
#include <ctldl/utility/contracts.hpp>

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
struct SparsityToFactorizeTridiagonal {
  SparsityViewStructural diag;
  SparsityViewStructural subdiag;
  PermutationViewStructural permutation;

  constexpr bool has_consistent_dim() const {
    return all_equal({diag.numRows(), diag.numCols(), subdiag.numRows(),
                      subdiag.numCols(), permutation.size()});
  }

  constexpr std::size_t dim() const {
    pre(has_consistent_dim());
    return diag.numRows();
  }
};

}  // namespace ctldl
