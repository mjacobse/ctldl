#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>

#include <array>

namespace ctldl {

template <class FactorData>
struct EmptyFactorDataLeft {
  static constexpr auto num_rows = FactorData::sparsity.num_rows;
  using Value = typename FactorData::Value;
  static constexpr auto sparsity = makeEmptySparsityCSR<num_rows, 0>();
  static constexpr std::array<Value, 0> L{};
  static constexpr auto permutation_row = FactorData::permutation_row;
  static constexpr Permutation<0> permutation_col{};
};

}  // namespace ctldl
