#pragma once

#include <ctldl/sparsity/sparsity.hpp>

namespace ctldl {

template <class T>
struct LowerTriangleBlocked3x3 {
  T block11;
  T block21;
  T block22;
  T block31;
  T block32;
  T block33;
};

consteval auto defineStaticSparsity(
    const LowerTriangleBlocked3x3<SparsityDynamic>& sparsity) {
  return LowerTriangleBlocked3x3{
      defineStaticSparsity(sparsity.block11),
      defineStaticSparsity(sparsity.block21),
      defineStaticSparsity(sparsity.block22),
      defineStaticSparsity(sparsity.block31),
      defineStaticSparsity(sparsity.block32),
      defineStaticSparsity(sparsity.block33),
  };
}

template <class T>
struct LowerTriangleBlocked {
  T block11;
  T block21;
  T block22;
};

consteval auto defineStaticSparsity(
    const LowerTriangleBlocked<SparsityDynamic>& sparsity) {
  return LowerTriangleBlocked{
      defineStaticSparsity(sparsity.block11),
      defineStaticSparsity(sparsity.block21),
      defineStaticSparsity(sparsity.block22),
  };
}

}  // namespace ctldl
