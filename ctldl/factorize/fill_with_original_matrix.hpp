#pragma once

#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t entry_index, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValue(
    const Matrix& input, FactorData& fact) {
  using Value = typename FactorData::Value;
  constexpr auto entry_orig =
      FactorData::origEntry(FactorData::sparsity.entries()[entry_index]);
  fact.L[entry_index] = static_cast<Value>(
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input));
}

template <std::size_t i, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillDiagonalWithOriginalMatrixValue(
    const Matrix& input, FactorData& fact) {
  using Value = typename FactorData::Value;
  constexpr auto i_orig = FactorData::origRowIndex(i);
  fact.D[i] = static_cast<Value>(getMatrixValueAt<i_orig, i_orig>(input));
}

template <std::size_t... EntryIndices, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValues(
    const Matrix& input, FactorData& fact,
    std::index_sequence<EntryIndices...>) {
  (fillRowWithOriginalMatrixValue<EntryIndices>(input, fact), ...);
}

template <std::size_t i, class Matrix, class FactorData>
[[gnu::always_inline]] inline auto fillRowWithOriginalMatrixValues(
    const Matrix& input, FactorData& fact) {
  constexpr auto row_begin = FactorData::sparsity.rowBeginIndices()[i];
  constexpr auto row_end = FactorData::sparsity.rowBeginIndices()[i + 1];
  fillRowWithOriginalMatrixValues(input, fact,
                                  makeIndexSequence<row_begin, row_end>());
}

template <std::size_t... RowIndices, class Matrix, class FactorData>
void fillWithOriginalMatrixValues(const Matrix& input, FactorData& fact,
                                  std::index_sequence<RowIndices...>) {
  (fillRowWithOriginalMatrixValues<RowIndices>(input, fact), ...);
}

template <std::size_t... RowIndices, class Matrix, class FactorData>
void fillWithOriginalMatrixValuesIncludingDiagonal(
    const Matrix& input, FactorData& fact, std::index_sequence<RowIndices...>) {
  (fillRowWithOriginalMatrixValues<RowIndices>(input, fact), ...);
  (fillDiagonalWithOriginalMatrixValue<RowIndices>(input, fact), ...);
}

template <class Matrix, class FactorData>
void fillWithOriginalMatrixValues(const Matrix& input, FactorData& fact) {
  constexpr auto num_rows = std::size_t{FactorData::sparsity.numRows()};
  fillWithOriginalMatrixValues(input, fact,
                               std::make_index_sequence<num_rows>());
}

template <class Matrix, class FactorData>
void fillWithOriginalMatrixValuesIncludingDiagonal(const Matrix& input,
                                                   FactorData& fact) {
  constexpr auto num_rows = std::size_t{FactorData::sparsity.numRows()};
  fillWithOriginalMatrixValuesIncludingDiagonal(
      input, fact, std::make_index_sequence<num_rows>());
}

}  // namespace ctldl
