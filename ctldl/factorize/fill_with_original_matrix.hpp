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
      FactorData::origEntry(FactorData::sparsity.entries[entry_index]);
  fact.L[entry_index] = static_cast<Value>(
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input));
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
  constexpr auto row_begin = FactorData::sparsity.row_begin_indices[i];
  constexpr auto row_end = FactorData::sparsity.row_begin_indices[i + 1];
  fillRowWithOriginalMatrixValues(input, fact,
                                  makeIndexSequence<row_begin, row_end>());
}

}  // namespace ctldl
