#pragma once

#include <ctldl/factorize/fill_with_original_matrix.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <cstddef>
#include <type_traits>

namespace ctldl {

template <std::size_t num_rows, std::size_t num_cols>
struct FactorInitNone {
  static constexpr auto sparsity =
      makeEmptySparsityStatic<num_rows, num_cols>();
};

template <std::size_t entry_index, class Init, class FactorData>
auto getInitialFactorValuePlainL(const Init& init, const FactorData& factor) {
  constexpr auto& sparsity = FactorData::sparsity;
  if constexpr (std::is_same_v<Init, FactorInitNone<sparsity.numRows(),
                                                    sparsity.numCols()>>) {
    return factor.L[entry_index];
  } else {
    constexpr auto i = std::size_t{sparsity.entries()[entry_index].row_index};
    constexpr auto j = std::size_t{sparsity.entries()[entry_index].col_index};
    constexpr auto entry_orig = FactorData::origEntry({i, j});
    return getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(init);
  }
}

template <std::size_t entry_index, class Init, class FactorData>
auto getInitialFactorValueL(const Init& init, const FactorData& factor) {
  using Value = typename FactorData::Value;
  return static_cast<Value>(
      getInitialFactorValuePlainL<entry_index>(init, factor));
}

template <std::size_t row_index, class Init, class FactorData>
auto getInitialFactorValuePlainD(const Init& init, const FactorData& factor) {
  constexpr auto& sparsity = FactorData::sparsity;
  if constexpr (std::is_same_v<Init, FactorInitNone<sparsity.numRows(),
                                                    sparsity.numCols()>>) {
    return factor.D[row_index];
  } else {
    constexpr auto row_index_orig = FactorData::origRowIndex(row_index);
    return getMatrixValueAt<row_index_orig, row_index_orig>(init);
  }
}

template <std::size_t row_index, class Init, class FactorData>
auto getInitialFactorValueD(const Init& init, const FactorData& factor) {
  using Value = typename FactorData::Value;
  return static_cast<Value>(
      getInitialFactorValuePlainD<row_index>(init, factor));
}

template <std::size_t row_index, class Init, class FactorData>
void initializeFactorRow(const Init& init, FactorData& factor) {
  if constexpr (!std::is_same_v<
                    Init, FactorInitNone<FactorData::sparsity.numRows(),
                                         FactorData::sparsity.numCols()>>) {
    fillRowWithOriginalMatrixValues<row_index>(init, factor);
  }
}

}  // namespace ctldl
