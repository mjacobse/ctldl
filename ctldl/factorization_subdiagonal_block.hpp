#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/utility/square.hpp>

#include <array>
#include <cstddef>

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

  template <class Matrix, class DiagonalBlock>
  void factor(const Matrix& matrix, const DiagonalBlock& above) {
    factorImpl(matrix, above);
  }

  std::array<Value, nnz> L;

 private:
  template <std::size_t entry_index_ij = 0, class Matrix, class DiagonalBlock>
  [[gnu::always_inline]] void factorImpl(const Matrix& matrix,
                                         const DiagonalBlock& above) {
    using DiagonalBlockSparsity = typename DiagonalBlock::Sparsity;

    if constexpr (entry_index_ij < nnz) {
      constexpr auto i = std::size_t{Sparsity::entries[entry_index_ij].row_index};
      constexpr auto j = std::size_t{Sparsity::entries[entry_index_ij].col_index};

      constexpr auto i_orig = permutation_row[i];
      constexpr auto j_orig = permutation_col[j];
      auto Lij = getMatrixValueAt<i_orig, j_orig>(matrix);

      static constexpr auto contributions =
          getContributionsMixed<DiagonalBlockSparsity, Sparsity, i, j>();
      for (const auto c : contributions) {
        Lij -= L[c.entry_index_ik] * above.L[c.entry_index_jk] * above.D[c.k];
      }
      L[entry_index_ij] = Lij / above.D[j];

      factorImpl<entry_index_ij + 1>(matrix, above);
    }
  }
};

}  // namespace ctldl
