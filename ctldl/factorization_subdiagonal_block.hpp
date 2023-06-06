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

  template <class Matrix, class DiagonalBlock>
  void factor(const Matrix& matrix, const DiagonalBlock& above) {
    factorImpl(matrix, above);
  }

  std::array<Value, nnz> L;

 private:
  template <std::size_t entry_index_ij, class Matrix, class DiagonalBlock>
  [[gnu::always_inline]] void factorImplEntry(const Matrix& matrix,
                                              const DiagonalBlock& above) {
    using DiagonalBlockSparsity = typename DiagonalBlock::Sparsity;

    constexpr auto i = std::size_t{Sparsity::entries[entry_index_ij].row_index};
    constexpr auto j = std::size_t{Sparsity::entries[entry_index_ij].col_index};

    constexpr auto entry_orig =
        permutedEntry(Entry{i, j}, permutation_row, permutation_col);
    auto Lij =
        getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(matrix);

    static constexpr auto contributions =
        getContributionsMixed<DiagonalBlockSparsity, Sparsity, i, j>();
    for (const auto c : contributions) {
      Lij -= L[c.entry_index_ik] * above.L[c.entry_index_jk] * above.D[c.k];
    }
    L[entry_index_ij] = Lij / above.D[j];
  }

  template <std::size_t... EntryIndices, class Matrix, class DiagonalBlock>
  [[gnu::always_inline]] void factorImplRow(
      const Matrix& matrix, const DiagonalBlock& above,
      std::index_sequence<EntryIndices...>) {
    (factorImplEntry<EntryIndices>(matrix, above), ...);
  }

  template <std::size_t i, class Matrix, class DiagonalBlock>
  [[gnu::always_inline]] void factorImplRow(const Matrix& matrix,
                                            const DiagonalBlock& above) {
    constexpr auto row_begin = std::size_t{Sparsity::row_begin_indices[i]};
    constexpr auto row_end = std::size_t{Sparsity::row_begin_indices[i + 1]};
    factorImplRow(matrix, above, makeIndexSequence<row_begin, row_end>());
  }

  template <std::size_t... RowIndices, class Matrix, class DiagonalBlock>
  void factorImpl(const Matrix& matrix, const DiagonalBlock& above,
                  std::index_sequence<RowIndices...>) {
    (factorImplRow<RowIndices>(matrix, above), ...);
  }

  template <class Matrix, class DiagonalBlock>
  void factorImpl(const Matrix& matrix, const DiagonalBlock& above) {
    factorImpl(matrix, above, std::make_index_sequence<num_rows>());
  }
};

}  // namespace ctldl
