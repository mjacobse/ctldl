#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>

#include <array>

namespace ctldl {

template <class FactorData>
struct EmptyFactorDataDiagonal {
  using Value = typename FactorData::Value;
  static constexpr std::array<Value, 0> D{};
};

}  // namespace ctldl
