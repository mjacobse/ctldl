#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/symbolic/filled_in_sparsity_repeating_arrowhead.hpp>
#include <ctldl/symbolic/repeating_block_tridiagonal.hpp>

#include <cstddef>

namespace ctldl {

template <auto sparsity_in_A, auto sparsity_in_B, auto permutation>
constexpr auto getFilledInSparsityRepeating() {
  static_assert(isSquare(sparsity_in_A));
  static_assert(isSquare(sparsity_in_B));
  static_assert(sparsity_in_B.numRows() == sparsity_in_A.numRows());
  constexpr auto dim = std::size_t{sparsity_in_A.numRows()};

  const auto sparsity_filled = getFilledInSparsityRepeatingArrowhead<
      sparsity_in_A, sparsity_in_B, makeEmptySparsity<0, dim>(), permutation,
      PermutationStatic<0>{}>();
  return RepeatingBlockTridiagonal{sparsity_filled.diag,
                                   sparsity_filled.subdiag};
};

}  // namespace ctldl
