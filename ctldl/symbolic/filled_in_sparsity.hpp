#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/symbolic/get_entries_with_fill.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>

namespace ctldl {

constexpr auto getFilledInSparsity(const SparsityView sparsity,
                                   const PermutationView permutation) {
  pre(isSquare(sparsity));
  pre(permutation.size() == sparsity.numRows());
  const auto dim = std::size_t{sparsity.numRows()};
  return SparsityDynamic(
      dim, dim,
      getEntriesWithFill(SparsityDynamicCSR(
          getSparsityDynamicLowerTriangle(sparsity, permutation))));
}

}  // namespace ctldl
