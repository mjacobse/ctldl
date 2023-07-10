#pragma once

#include <ctldl/empty_factor_data_diagonal.hpp>
#include <ctldl/empty_factor_data_left.hpp>
#include <ctldl/empty_matrix_input.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/permutation/permuted_entry_lower_triangle.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/utility/make_index_sequence.hpp>
#include <ctldl/utility/square.hpp>

#include <array>
#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t entry_index_ij, class FactorData, class Matrix,
          class FactorDataAbove>
[[gnu::always_inline]] inline void factorEntryWiseSubdiagonalImplRow(
    FactorData& self, const Matrix& matrix, const FactorDataAbove& above) {
  constexpr auto& sparsity = FactorData::sparsity;
  constexpr auto& sparsity_above = FactorDataAbove::sparsity;
  using Value = typename FactorData::Value;

  constexpr auto i = std::size_t{sparsity.entries[entry_index_ij].row_index};
  constexpr auto j = std::size_t{sparsity.entries[entry_index_ij].col_index};

  constexpr auto entry_orig = permutedEntry(
      Entry{i, j}, FactorData::permutation_row, FactorData::permutation_col);
  auto Lij = static_cast<Value>(
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(matrix));

  static constexpr auto contributions =
      getContributionsMixed<sparsity_above, sparsity, i, j>();
  for (const auto c : contributions) {
    Lij -= self.L[c.entry_index_ik] * above.L[c.entry_index_jk] * above.D[c.k];
  }
  self.L[entry_index_ij] = Lij / above.D[j];
}

template <std::size_t... EntryIndices, class FactorData, class Matrix,
          class FactorDataAbove>
[[gnu::always_inline]] inline void factorEntryWiseSubdiagonalImplRow(
    FactorData& self, const Matrix& matrix, const FactorDataAbove& above,
    std::index_sequence<EntryIndices...>) {
  (factorEntryWiseSubdiagonalImplRow<EntryIndices>(self, matrix, above), ...);
}

template <std::size_t i, class FactorData, class Matrix, class FactorDataAbove>
[[gnu::always_inline]] inline void factorizeEntryWiseSubdiagonalImpl(
    FactorData& self, const Matrix& matrix, const FactorDataAbove& above) {
  constexpr auto& sparsity = FactorData::sparsity;
  constexpr auto row_begin = std::size_t{sparsity.row_begin_indices[i]};
  constexpr auto row_end = std::size_t{sparsity.row_begin_indices[i + 1]};
  factorEntryWiseSubdiagonalImplRow(self, matrix, above,
                                    makeIndexSequence<row_begin, row_end>());
}

template <std::size_t... RowIndices, class FactorData, class Matrix,
          class FactorDataAbove>
void factorizeEntryWiseSubdiagonalImpl(FactorData& self,
                                       const Matrix& input,
                                       const FactorDataAbove& above,
                                       std::index_sequence<RowIndices...>) {
  (factorizeEntryWiseSubdiagonalImpl<RowIndices>(self, input, above), ...);
}

template <std::size_t i, class FactorData, class FactorDataDiag>
[[gnu::always_inline]] inline auto applyContributionsRowDiagonal(
    const FactorData& fact, const FactorDataDiag& diag,
    const typename FactorData::Value value_init) {
  constexpr auto& sparsity = FactorData::sparsity;

  constexpr auto row_begin = std::size_t{sparsity.row_begin_indices[i]};
  constexpr auto row_end = std::size_t{sparsity.row_begin_indices[i + 1]};
  auto value = value_init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = sparsity.entries[entry_index_ij].col_index;
    value -= square(fact.L[entry_index_ij]) * diag.D[j];
  }
  return value;
}

template <class FactorData, class FactorDataDiag, std::size_t num_contributions>
[[gnu::always_inline]] inline auto applyContributions(
    const FactorData& fact, const FactorDataDiag& diag,
    const std::array<Contribution, num_contributions>& contributions,
    const typename FactorData::Value value_init) {
  auto value = value_init;
  for (const auto c : contributions) {
    value -= fact.L[c.entry_index_ik] * fact.L[c.entry_index_jk] * diag.D[c.k];
  }
  return value;
}

template <std::size_t entry_index_ij, class FactorData, class Matrix,
          class FactorDataLeft, class FactorDataAbove>
[[gnu::always_inline]] inline auto factorizeEntryWiseImplRow(
    FactorData& self, const Matrix& input, const FactorDataLeft& left,
    const FactorDataAbove& above) {
  constexpr auto& sparsity = FactorData::sparsity;
  constexpr auto& sparsity_left = FactorDataLeft::sparsity;
  using Value = typename FactorData::Value;

  constexpr auto i = std::size_t{sparsity.entries[entry_index_ij].row_index};
  constexpr auto j = std::size_t{sparsity.entries[entry_index_ij].col_index};

  constexpr auto entry_orig =
      permutedEntryLowerTriangle(Entry{i, j}, FactorData::permutation);
  auto Lij = static_cast<Value>(
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input));

  static constexpr auto contributions_left =
      getContributionsRectangular<sparsity_left, i, j>();
  static constexpr auto contributions_self =
      getContributionsLowerTriangle<sparsity, i, j>();

  Lij = applyContributions(left, above, contributions_left, Lij);
  Lij = applyContributions(self, self, contributions_self, Lij);
  Lij /= self.D[j];
  self.L[entry_index_ij] = Lij;

  return square(Lij) * self.D[j];
}

template <std::size_t... EntryIndices, class FactorData, class Matrix,
          class FactorDataLeft, class FactorDataAbove>
[[gnu::always_inline]] inline auto factorizeEntryWiseImplRow(
    FactorData& self, const Matrix& input, const FactorDataLeft& left,
    const FactorDataAbove& above, const typename FactorData::Value Di_init,
    std::index_sequence<EntryIndices...>) {
  return (Di_init - ... -
          factorizeEntryWiseImplRow<EntryIndices>(self, input, left, above));
}

template <std::size_t i, class FactorData, class Matrix, class FactorDataLeft,
          class FactorDataAbove>
[[gnu::always_inline]] inline void factorizeEntryWiseImpl(
    FactorData& self, const Matrix& input, const FactorDataLeft& left,
    const FactorDataAbove& above) {
  constexpr auto& sparsity = FactorData::sparsity;
  constexpr auto& sparsity_left = FactorDataLeft::sparsity;
  static_assert(sparsity.num_rows == sparsity_left.num_rows);
  using Value = typename FactorData::Value;

  constexpr auto row_begin = std::size_t{sparsity.row_begin_indices[i]};
  constexpr auto row_end = std::size_t{sparsity.row_begin_indices[i + 1]};

  constexpr auto i_orig = FactorData::permutation[i];
  auto Di = static_cast<Value>(getMatrixValueAt<i_orig, i_orig>(input));
  Di = applyContributionsRowDiagonal<i>(left, above, Di);
  Di = factorizeEntryWiseImplRow(self, input, left, above, Di,
                                 makeIndexSequence<row_begin, row_end>());
  self.D[i] = Di;
}

template <std::size_t... RowIndices, class FactorData, class Matrix,
          class FactorDataLeft, class FactorDataAbove>
void factorizeEntryWiseImpl(FactorData& self, const Matrix& input,
                            const FactorDataLeft& left,
                            const FactorDataAbove& above,
                            std::index_sequence<RowIndices...>) {
  (factorizeEntryWiseImpl<RowIndices>(self, input, left, above), ...);
}

template <class FactorDataAbove, class MatrixLeft, class MatrixSelf,
          class FactorDataLeft, class FactorData>
void factorizeEntryWise(const FactorDataAbove& above,
                        const MatrixLeft& input_left,
                        const MatrixSelf& input_self, FactorDataLeft& left,
                        FactorData& self) {
  constexpr auto num_rows = std::size_t{FactorData::sparsity.num_rows};
  factorizeEntryWiseSubdiagonalImpl(left, input_left, above,
                                    std::make_index_sequence<num_rows>());
  factorizeEntryWiseImpl(self, input_self, left, above,
                         std::make_index_sequence<num_rows>());
}

template <class FactorData, class Matrix>
void factorizeEntryWise(FactorData& self, const Matrix& input) {
  constexpr EmptyFactorDataLeft<FactorData> empty_left;
  constexpr EmptyFactorDataDiagonal<decltype(empty_left)> empty_above;
  constexpr EmptyMatrixInput empty_input_left;
  factorizeEntryWise(empty_above, empty_input_left, input, empty_left, self);
}

template <class FactorData, class Matrix>
void factorize(FactorData& self, const Matrix& matrix,
               FactorizeMethodEntryWise) {
  factorizeEntryWise(self, matrix);
}

}  // namespace ctldl
