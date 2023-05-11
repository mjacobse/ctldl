#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/symbolic/get_entries_with_fill.hpp>

#include <cassert>

namespace ctldl {

template <class Sparsity, class PermutationIn = PermutationIdentity>
struct FilledInSparsity {
  static_assert(Sparsity::num_rows == Sparsity::num_cols);
  static constexpr auto num_rows = Sparsity::num_rows;
  static constexpr auto num_cols = Sparsity::num_cols;
  static constexpr auto dim = num_rows;
  static constexpr Permutation<dim> permutation{PermutationIn::permutation};
  static constexpr auto entries = getEntriesWithFill<
      SparsityCSR<SparsityLowerTriangle<Sparsity, PermutationIn>>>();
};

}  // namespace ctldl
