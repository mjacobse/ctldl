#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <std::size_t dim>
struct EmptyFactorSparsityDiagonal {
  static constexpr auto num_rows = dim;
  static constexpr auto num_cols = dim;
  static constexpr std::array<std::size_t, num_rows + 1> row_begin_indices{};
  static constexpr std::array<Entry, 0> entries{};
};

template <class FactorData>
struct EmptyFactorDataDiagonal {
  using Value = typename FactorData::Value;
  using Sparsity = EmptyFactorSparsityDiagonal<FactorData::Sparsity::num_cols>;
  static constexpr std::array<Value, 0> D{};
};

}  // namespace ctldl
