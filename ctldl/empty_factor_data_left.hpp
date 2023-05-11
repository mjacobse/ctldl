#pragma once

#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cassert>
#include <cstddef>

namespace ctldl {

template <std::size_t num_rows_>
struct EmptyFactorSparsityLeft {
  static constexpr auto num_rows = num_rows_;
  static constexpr auto num_cols = std::size_t{0};
  static constexpr std::array<std::size_t, num_rows + 1> row_begin_indices{};
  static constexpr std::array<Entry, 0> entries{};

  static constexpr bool isNonZero(const std::size_t /*i*/,
                                  const std::size_t /*j*/) {
    return false;
  }

  static constexpr std::size_t entryIndex(const std::size_t /*i*/,
                                          const std::size_t /*j*/) {
    assert(false);
    return 0;
  }
};

template <class FactorData>
struct EmptyFactorDataLeft {
  static constexpr auto num_rows = FactorData::Sparsity::num_rows;
  using Value = typename FactorData::Value;
  using Sparsity = EmptyFactorSparsityLeft<num_rows>;
  static constexpr std::array<Value, 0> L{};
  static constexpr std::array<Value, 0> D{};
  static constexpr auto permutation_row = FactorData::permutation_row;
  static constexpr Permutation<0> permutation_col{};
};

}  // namespace ctldl
