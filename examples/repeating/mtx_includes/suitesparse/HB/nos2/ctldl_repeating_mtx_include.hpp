#pragma once

#include <ctldl/factor_data/sparsity_to_factorize_link.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_outer.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_start.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <array>
#include <cstddef>

constexpr auto getRepeatingMtxSparsityStart() {
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsity<3, 0>();
  constexpr auto outer = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeStart{diag, next, outer, permutation};
}

constexpr auto getRepeatingMtxSparsityTridiag() {
  constexpr auto diag = ctldl::makeSparsity<3, 3>({{0, 0}, {1, 1}, {2, 2}});
  constexpr auto subdiag =
      ctldl::makeSparsity<3, 3>({{0, 0}, {1, 1}, {2, 1}, {1, 2}, {2, 2}});
  constexpr auto permutation = ctldl::Permutation<3>{{0, 1, 2}};
  return ctldl::SparsityToFactorizeTridiagonal{diag, subdiag, permutation};
}

constexpr auto getRepeatingMtxSparsityLink() {
  constexpr auto prev = ctldl::makeEmptySparsity<0, 3>();
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeLink{prev, diag, next, permutation};
}

constexpr auto getRepeatingMtxSparsityOuter() {
  constexpr auto subdiag = ctldl::makeEmptySparsity<0, 3>();
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeOuter{subdiag, diag, permutation};
}
