#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>

#include <array>

namespace ctldl {

template <class FactorData>
struct EmptyFactorDataDiagonal {
  static constexpr auto sparsity = makeEmptySparsityStaticCSR<0, 0>();
  using Value = typename FactorData::Value;
  static constexpr std::array<Value, 0> D{};
  static constexpr std::array<Value, 0> L{};
};

}  // namespace ctldl
