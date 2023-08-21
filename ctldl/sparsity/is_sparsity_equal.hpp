#pragma once

#include <ctldl/sparsity/is_sparsity_subset.hpp>

namespace ctldl {

template <class SparsityLhs, class SparsityRhs>
constexpr bool isSparsityEqual(
    const SparsityLhs& sparsity_lhs, const SparsityRhs& sparsity_rhs) {
  return isSparsitySubset(sparsity_lhs, sparsity_rhs) &&
         isSparsitySubset(sparsity_rhs, sparsity_lhs);
}

}  // namespace ctldl
