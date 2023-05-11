#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/get_entries.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/symbolic/get_filled_in_is_nonzero_info.hpp>

#include <cassert>

namespace ctldl {

template <class Sparsity, class PermutationIn = PermutationIdentity>
struct FilledInSparsity {
  static_assert(Sparsity::num_rows == Sparsity::num_cols);
  static constexpr auto num_rows = Sparsity::num_rows;
  static constexpr auto num_cols = Sparsity::num_cols;
  static constexpr auto dim = num_rows;
  static constexpr Permutation<dim> permutation{PermutationIn::permutation};
  static constexpr auto is_nonzero = getFilledInIsNonzeroInfo<
      SparsityCSR<SparsityLowerTriangle<Sparsity, PermutationIn>>>();
  static constexpr auto entries = getEntries([] { return is_nonzero; });
};

}  // namespace ctldl
