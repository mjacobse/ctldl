#pragma once

#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cassert>
#include <cstddef>

namespace ctldl {

template <std::size_t num_rows_>
struct EmptyFactorSparsity {
  static constexpr auto num_rows = num_rows_;
  static constexpr auto num_cols = std::size_t{0};
  static constexpr std::array<std::size_t, num_rows + 1> row_begin_indices{};
  static constexpr std::array<Entry, 0> entries{};
  static constexpr std::array<std::array<bool, 0>, num_rows> is_nonzero{};

  static constexpr std::size_t entryIndex(const std::size_t /*i*/,
                                          const std::size_t /*j*/) {
    assert(false);
    return 0;
  }
};

template <std::size_t num_rows_, class Value_>
struct EmptyFactorData {
  static constexpr auto num_rows = num_rows_;
  using Value = Value_;
  using Sparsity = EmptyFactorSparsity<num_rows>;
  static constexpr std::array<Value, 0> L{};
  static constexpr std::array<Value, 0> D{};
};

}  // namespace ctldl
