#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/symbolic/get_entries_with_fill.hpp>

#include <cstddef>

namespace ctldl {

template <auto sparsity, auto permutation = Permutation<sparsity.num_rows>{}>
constexpr auto getFilledInSparsity() {
  static_assert(sparsity.num_rows == sparsity.num_cols);
  constexpr auto dim = std::size_t{sparsity.num_rows};
  return makeSparsity<dim, dim>(
      getEntriesWithFill<SparsityCSR(
          getSparsityLowerTriangle<sparsity>(permutation))>());
}

}  // namespace ctldl
