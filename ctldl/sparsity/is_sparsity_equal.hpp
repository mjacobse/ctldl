#pragma once

#include <ctldl/sparsity/is_sparsity_subset.hpp>
#include <ctldl/sparsity/sparsity.hpp>

namespace ctldl {

constexpr bool isSparsityEqual(const SparsityView sparsity_lhs,
                               const SparsityView sparsity_rhs) {
  return isSparsitySubset(sparsity_lhs, sparsity_rhs) &&
         isSparsitySubset(sparsity_rhs, sparsity_lhs);
}

}  // namespace ctldl
