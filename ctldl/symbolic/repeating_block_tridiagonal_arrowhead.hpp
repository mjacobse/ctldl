#pragma once

#include <ctldl/sparsity/sparsity.hpp>

namespace ctldl {

template <class T>
struct RepeatingBlockTridiagonalArrowhead {
  T diag;
  T subdiag;
  T outer;
};

consteval auto defineStaticSparsity(
    const RepeatingBlockTridiagonalArrowhead<SparsityDynamic>& sparsity) {
  return RepeatingBlockTridiagonalArrowhead{
      defineStaticSparsity(sparsity.diag),
      defineStaticSparsity(sparsity.subdiag),
      defineStaticSparsity(sparsity.outer),
  };
}

}  // namespace ctldl
