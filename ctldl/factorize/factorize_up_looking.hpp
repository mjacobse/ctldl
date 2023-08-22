#pragma once

#include <ctldl/factor_data/empty_factor_data_diagonal.hpp>
#include <ctldl/factor_data/empty_factor_data_left.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/matrix/empty_matrix_input.hpp>
#include <ctldl/permutation/permuted_entry.hpp>
#include <ctldl/sparsity/get_influenced.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/sparsity/is_sparsity_subset.hpp>
#include <ctldl/symbolic/is_chordal_blocked.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t entry_index, class FactorData>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerSelf(
    FactorData& fact) {
  constexpr auto i = FactorData::sparsity.entries[entry_index].row_index;
  constexpr auto j = FactorData::sparsity.entries[entry_index].col_index;

  const auto Lij_scaled = fact.L[entry_index];

  static constexpr auto influenced_list =
      getInfluencedList<i, j, FactorData::sparsity, FactorData::sparsity>(
          [](const std::size_t i, const std::size_t /*j*/) { return i; });
  for (const auto influenced : influenced_list) {
    fact.L[influenced.entry_index_target] -=
        fact.L[influenced.entry_index_source] * Lij_scaled;
  }

  const auto Lij = Lij_scaled / fact.D[j];
  fact.L[entry_index] = Lij;
  return Lij_scaled * Lij;
}

template <std::size_t... EntryIndices, class FactorData>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerSelf(
    FactorData& fact, typename FactorData::Value Di_init,
    std::index_sequence<EntryIndices...>) {
  return (Di_init - ... - factorizeUpLookingInnerSelf<EntryIndices>(fact));
}

template <std::size_t entry_index, class FactorDataAbove, class FactorData,
          class FactorDataLeft>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorDataAbove& above, FactorDataLeft& left, FactorData& self) {
  constexpr auto i = FactorDataLeft::sparsity.entries[entry_index].row_index;
  constexpr auto j = FactorDataLeft::sparsity.entries[entry_index].col_index;

  const auto Lij_scaled = left.L[entry_index];

  static constexpr auto influenced_list_left =
      getInfluencedList<i, j, FactorData::sparsity, FactorDataLeft::sparsity>(
          [](const std::size_t /*i*/, const std::size_t /*j*/) {
            return FactorData::sparsity.num_rows;
          });
  for (const auto influenced : influenced_list_left) {
    left.L[influenced.entry_index_target] -=
        above.L[influenced.entry_index_source] * Lij_scaled;
  }
  static constexpr auto influenced_list_self =
      getInfluencedList<i, j, FactorDataLeft::sparsity, FactorData::sparsity>(
          [](const std::size_t i, const std::size_t /*j*/) { return i; });
  for (const auto influenced : influenced_list_self) {
    self.L[influenced.entry_index_target] -=
        left.L[influenced.entry_index_source] * Lij_scaled;
  }

  const auto Lij = Lij_scaled / above.D[j];
  left.L[entry_index] = Lij;
  return Lij_scaled * Lij;
}

template <std::size_t... EntryIndices, class FactorDataAbove,
          class FactorDataLeft, class FactorData>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorDataAbove& above, FactorDataLeft& left, FactorData& self,
    typename FactorData::Value Di_init, std::index_sequence<EntryIndices...>) {
  return (Di_init - ... -
          factorizeUpLookingInnerLeft<EntryIndices>(above, left, self));
}

template <std::size_t entry_index, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValueSubdiagonal(
    const Matrix& input, FactorData& fact) {
  using Value = typename FactorData::Value;
  constexpr auto entry_orig =
      permutedEntry(FactorData::sparsity.entries[entry_index],
                    FactorData::permutation_row, FactorData::permutation_col);
  fact.L[entry_index] = static_cast<Value>(
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input));
}

template <std::size_t... EntryIndices, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValuesSubdiagonal(
    const Matrix& input, FactorData& fact,
    std::index_sequence<EntryIndices...>) {
  (fillRowWithOriginalMatrixValueSubdiagonal<EntryIndices>(input, fact), ...);
}

template <std::size_t i, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValuesSubdiagonal(
    const Matrix& input, FactorData& fact) {
  constexpr auto row_begin = FactorData::sparsity.row_begin_indices[i];
  constexpr auto row_end = FactorData::sparsity.row_begin_indices[i + 1];
  fillRowWithOriginalMatrixValuesSubdiagonal(
      input, fact, makeIndexSequence<row_begin, row_end>());
}

template <std::size_t entry_index, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValueDiagonal(
    const Matrix& input, FactorData& fact) {
  using Value = typename FactorData::Value;
  constexpr auto entry_orig = permutedEntryLowerTriangle(
      FactorData::sparsity.entries[entry_index], FactorData::permutation);
  fact.L[entry_index] = static_cast<Value>(
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input));
}

template <std::size_t... EntryIndices, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValuesDiagonal(
    const Matrix& input, FactorData& fact,
    std::index_sequence<EntryIndices...>) {
  (fillRowWithOriginalMatrixValueDiagonal<EntryIndices>(input, fact), ...);
}

template <std::size_t i, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValuesDiagonal(
    const Matrix& input, FactorData& fact) {
  constexpr auto row_begin = FactorData::sparsity.row_begin_indices[i];
  constexpr auto row_end = FactorData::sparsity.row_begin_indices[i + 1];
  fillRowWithOriginalMatrixValuesDiagonal(
      input, fact, makeIndexSequence<row_begin, row_end>());
}

template <std::size_t i, class FactorDataAbove, class MatrixLeft,
          class MatrixSelf, class FactorDataLeft, class FactorData>
[[gnu::always_inline]] inline void factorizeUpLookingImpl(
    const FactorDataAbove& above, const MatrixLeft& input_left,
    const MatrixSelf& input_self, FactorDataLeft& left, FactorData& self) {
  constexpr auto& sparsity = FactorData::sparsity;
  constexpr auto& sparsity_left = FactorDataLeft::sparsity;
  static_assert(sparsity.num_rows == sparsity_left.num_rows);
  using Value = typename FactorData::Value;

  constexpr auto i_orig = FactorData::permutation[i];
  auto Di = static_cast<Value>(getMatrixValueAt<i_orig, i_orig>(input_self));

  fillRowWithOriginalMatrixValuesSubdiagonal<i>(input_left, left);
  fillRowWithOriginalMatrixValuesDiagonal<i>(input_self, self);

  constexpr auto row_begin_left = sparsity_left.row_begin_indices[i];
  constexpr auto row_end_left = sparsity_left.row_begin_indices[i + 1];
  Di = factorizeUpLookingInnerLeft(
      above, left, self, Di, makeIndexSequence<row_begin_left, row_end_left>());
  constexpr auto row_begin_self = sparsity.row_begin_indices[i];
  constexpr auto row_end_self = sparsity.row_begin_indices[i + 1];
  Di = factorizeUpLookingInnerSelf(
      self, Di, makeIndexSequence<row_begin_self, row_end_self>());
  self.D[i] = Di;
}

template <std::size_t... RowIndices, class FactorDataAbove, class MatrixLeft,
          class MatrixSelf, class FactorDataLeft, class FactorData>
void factorizeUpLookingImpl(const FactorDataAbove& above,
                            const MatrixLeft& input_left,
                            const MatrixSelf& input_self, FactorDataLeft& left,
                            FactorData& self,
                            std::index_sequence<RowIndices...>) {
  (factorizeUpLookingImpl<RowIndices>(above, input_left, input_self, left,
                                      self),
   ...);
}

template <class FactorDataAbove, class MatrixLeft, class MatrixSelf,
          class FactorDataLeft, class FactorData>
void factorizeUpLooking(const FactorDataAbove& above,
                        const MatrixLeft& input_left,
                        const MatrixSelf& input_self, FactorDataLeft& left,
                        FactorData& self) {
  static_assert(isChordalBlocked(FactorDataAbove::sparsity,
                                 FactorDataLeft::sparsity,
                                 FactorData::sparsity));
  static_assert(isSparsitySubsetLowerTriangle<MatrixSelf::sparsity>(
      FactorData::sparsity, FactorData::permutation));
  static_assert(isSparsitySubset(MatrixLeft::sparsity, FactorDataLeft::sparsity,
                                 FactorDataLeft::permutation_row,
                                 FactorDataLeft::permutation_col));
  constexpr auto num_rows = std::size_t{FactorData::sparsity.num_rows};
  factorizeUpLookingImpl(above, input_left, input_self, left, self,
                         std::make_index_sequence<num_rows>());
}

template <class FactorData, class Matrix>
void factorizeUpLooking(FactorData& self, const Matrix& matrix) {
  constexpr EmptyFactorDataLeft<FactorData> empty_left;
  constexpr EmptyFactorDataDiagonal<decltype(empty_left)> empty_above;
  constexpr EmptyMatrixInput<FactorData::sparsity.num_rows, 0> empty_input_left;
  factorizeUpLooking(empty_above, empty_input_left, matrix, empty_left, self);
}

template <class FactorData, class Matrix>
void factorize(FactorData& self, const Matrix& matrix,
               FactorizeMethodUpLooking) {
  factorizeUpLooking(self, matrix);
}

}  // namespace ctldl
