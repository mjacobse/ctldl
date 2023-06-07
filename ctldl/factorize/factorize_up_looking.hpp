#pragma once

#include <ctldl/permutation/permuted_entry.hpp>
#include <ctldl/sparsity/get_influenced.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t entry_index, class FactorData>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerSelf(
    FactorData& fact) {
  using Sparsity = FactorData::Sparsity;

  constexpr auto i = Sparsity::entries[entry_index].row_index;
  constexpr auto j = Sparsity::entries[entry_index].col_index;

  const auto Lij_scaled = fact.L[entry_index];

  static constexpr auto influenced_list =
      getInfluencedList<i, j, Sparsity, Sparsity>(
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

template <std::size_t entry_index, class FactorData, class FactorDataLeft>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorData& above, FactorDataLeft& left, FactorData& self) {
  using Sparsity = FactorData::Sparsity;
  using SparsityLeft = FactorDataLeft::Sparsity;

  constexpr auto i = SparsityLeft::entries[entry_index].row_index;
  constexpr auto j = SparsityLeft::entries[entry_index].col_index;

  const auto Lij_scaled = left.L[entry_index];

  static constexpr auto influenced_list_left =
      getInfluencedList<i, j, Sparsity, SparsityLeft>(
          [](const std::size_t /*i*/, const std::size_t /*j*/) {
            return Sparsity::num_rows;
          });
  for (const auto influenced : influenced_list_left) {
    left.L[influenced.entry_index_target] -=
        above.L[influenced.entry_index_source] * Lij_scaled;
  }
  static constexpr auto influenced_list_self =
      getInfluencedList<i, j, SparsityLeft, Sparsity>(
          [](const std::size_t i, const std::size_t /*j*/) { return i; });
  for (const auto influenced : influenced_list_self) {
    self.L[influenced.entry_index_target] -=
        left.L[influenced.entry_index_source] * Lij_scaled;
  }

  const auto Lij = Lij_scaled / above.D[j];
  left.L[entry_index] = Lij;
  return Lij_scaled * Lij;
}

template <std::size_t... EntryIndices, class FactorData, class FactorDataLeft>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorData& above, FactorDataLeft& left, FactorData& self,
    typename FactorData::Value Di_init, std::index_sequence<EntryIndices...>) {
  return (Di_init - ... -
          factorizeUpLookingInnerLeft<EntryIndices>(above, left, self));
}

template <std::size_t entry_index, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValueSubdiagonal(
    const Matrix& input, FactorData& fact) {
  using Sparsity = typename FactorData::Sparsity;
  constexpr auto entry_orig =
      permutedEntry(Sparsity::entries[entry_index], FactorData::permutation_row,
                    FactorData::permutation_col);
  fact.L[entry_index] =
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input);
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
  using Sparsity = typename FactorData::Sparsity;
  constexpr auto row_begin = Sparsity::row_begin_indices[i];
  constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
  fillRowWithOriginalMatrixValuesSubdiagonal(
      input, fact, makeIndexSequence<row_begin, row_end>());
}

template <std::size_t entry_index, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValueDiagonal(
    const Matrix& input, FactorData& fact) {
  using Sparsity = typename FactorData::Sparsity;
  constexpr auto entry_orig = permutedEntryLowerTriangle(
      Sparsity::entries[entry_index], FactorData::permutation);
  fact.L[entry_index] =
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input);
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
  using Sparsity = typename FactorData::Sparsity;
  constexpr auto row_begin = Sparsity::row_begin_indices[i];
  constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
  fillRowWithOriginalMatrixValuesDiagonal(
      input, fact, makeIndexSequence<row_begin, row_end>());
}

template <std::size_t i, class FactorData, class MatrixLeft, class MatrixSelf,
          class FactorDataLeft>
[[gnu::always_inline]] inline void factorizeUpLookingImpl(
    const FactorData& above, const MatrixLeft& input_left,
    const MatrixSelf& input_self, FactorDataLeft& left, FactorData& self) {
  using Sparsity = typename FactorData::Sparsity;
  using SparsityLeft = typename FactorDataLeft::Sparsity;
  static_assert(Sparsity::num_rows == SparsityLeft::num_rows);

  constexpr auto i_orig = FactorData::permutation[i];
  auto Di = getMatrixValueAt<i_orig, i_orig>(input_self);

  fillRowWithOriginalMatrixValuesSubdiagonal<i>(input_left, left);
  fillRowWithOriginalMatrixValuesDiagonal<i>(input_self, self);

  constexpr auto row_begin_left = SparsityLeft::row_begin_indices[i];
  constexpr auto row_end_left = SparsityLeft::row_begin_indices[i + 1];
  Di = factorizeUpLookingInnerLeft(
      above, left, self, Di, makeIndexSequence<row_begin_left, row_end_left>());
  constexpr auto row_begin_self = Sparsity::row_begin_indices[i];
  constexpr auto row_end_self = Sparsity::row_begin_indices[i + 1];
  Di = factorizeUpLookingInnerSelf(
      self, Di, makeIndexSequence<row_begin_self, row_end_self>());
  self.D[i] = Di;
}

template <std::size_t... RowIndices, class FactorData, class MatrixLeft,
          class MatrixSelf, class FactorDataLeft>
void factorizeUpLookingImpl(const FactorData& above,
                            const MatrixLeft& input_left,
                            const MatrixSelf& input_self, FactorDataLeft& left,
                            FactorData& self,
                            std::index_sequence<RowIndices...>) {
  (factorizeUpLookingImpl<RowIndices>(above, input_left, input_self, left,
                                      self),
   ...);
}

template <class FactorData, class MatrixLeft, class MatrixSelf,
          class FactorDataLeft>
void factorizeUpLooking(const FactorData& above, const MatrixLeft& input_left,
                        const MatrixSelf& input_self, FactorDataLeft& left,
                        FactorData& self) {
  constexpr auto num_rows = std::size_t{FactorData::Sparsity::num_rows};
  factorizeUpLookingImpl(above, input_left, input_self, left, self,
                         std::make_index_sequence<num_rows>());
}

}  // namespace ctldl
