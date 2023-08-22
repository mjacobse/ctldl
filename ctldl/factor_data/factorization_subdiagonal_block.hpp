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

template <auto sparsity_in, class Value_,
          auto permutation_row_in = PermutationIdentity{},
          auto permutation_col_in = PermutationIdentity{}>
class FactorizationSubdiagonalBlock {
 public:
  static constexpr auto sparsity = makeSparsityCSR(sparsity_in);
  static constexpr auto nnz = std::size_t{sparsity.nnz};
  static constexpr auto num_rows = std::size_t{sparsity.num_rows};
  static constexpr auto num_cols = std::size_t{sparsity.num_cols};
  using Value = Value_;
  static constexpr Permutation<num_rows> permutation_row{permutation_row_in};
  static constexpr Permutation<num_cols> permutation_col{permutation_col_in};

  std::array<Value, nnz> L;
};

}  // namespace ctldl
