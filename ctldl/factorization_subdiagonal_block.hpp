#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/permutation/permuted_entry.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

#include <array>
#include <cstddef>
#include <utility>

namespace ctldl {

template <class OriginalSparsity, class Value_,
          class PermutationRow = PermutationIdentity,
          class PermutationCol = PermutationIdentity>
class FactorizationSubdiagonalBlock {
 public:
  using Sparsity = SparsityCSR<OriginalSparsity>;
  static constexpr auto nnz = std::size_t{Sparsity::nnz};
  static constexpr auto num_rows = std::size_t{Sparsity::num_rows};
  static constexpr auto num_cols = std::size_t{Sparsity::num_cols};
  using Value = Value_;
  static constexpr Permutation<num_rows> permutation_row{
      PermutationRow::permutation};
  static constexpr Permutation<num_cols> permutation_col{
      PermutationCol::permutation};

  std::array<Value, nnz> L;
};

}  // namespace ctldl
