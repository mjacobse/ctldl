#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>

#include <array>

namespace ctldl {

template <class FactorData>
struct EmptyFactorDataLeft {
  static constexpr auto num_rows = FactorData::sparsity.numRows();
  using Value = typename FactorData::Value;
  static constexpr auto sparsity = makeEmptySparsityCSR<num_rows, 0>();
  static constexpr std::array<Value, 0> L{};
  static constexpr auto permutation_row = FactorData::permutation_row;
  static constexpr Permutation<0> permutation_col{};

  static constexpr auto origRowIndex(const std::size_t factor_row_index) {
    return std::size_t{permutation_row[factor_row_index]};
  }
  static constexpr auto origColIndex(const std::size_t factor_col_index) {
    return std::size_t{permutation_col[factor_col_index]};
  }
};

}  // namespace ctldl
