#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/permutation/permuted_entry.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

#include <array>
#include <cstddef>
#include <utility>

namespace ctldl {

template <Sparsity sparsity_in, class Value_,
          Permutation<sparsity_in.num_rows> permutation_row_in =
              PermutationIdentity{},
          Permutation<sparsity_in.num_cols> permutation_col_in =
              PermutationIdentity{}>
class FactorizationSubdiagonalBlock {
 public:
  static constexpr auto sparsity = makeSparsityCSR(sparsity_in);
  static constexpr auto nnz = std::size_t{sparsity.nnz};
  static constexpr auto num_rows = std::size_t{sparsity.num_rows};
  static constexpr auto num_cols = std::size_t{sparsity.num_cols};
  using Value = Value_;
  static constexpr auto permutation_row = permutation_row_in;
  static constexpr auto permutation_col = permutation_col_in;

  /**
   * Given an entry location in the factor, returns the location of the
   * corresponding entry in the original matrix, i.e. the location before the
   * permutation for the factorization was applied.
   */
  static constexpr Entry origEntry(const Entry factor_entry) {
    return permutedEntry(factor_entry, permutation_row, permutation_col);
  };
  /**
   * For a given row index in the factor, returns the corresponding row index in
   * the original matrix, i.e. the row index before the permutation for the
   * factorization was applied.
   */
  static constexpr auto origRowIndex(const std::size_t factor_row_index) {
    return std::size_t{permutation_row[factor_row_index]};
  }
  /**
   * For a given column index in the factor, returns the corresponding column
   * index in the original matrix, i.e. the column index before the permutation
   * for the factorization was applied.
   */
  static constexpr auto origColIndex(const std::size_t factor_col_index) {
    return std::size_t{permutation_col[factor_col_index]};
  }

  std::array<Value, nnz> L;
};

}  // namespace ctldl
