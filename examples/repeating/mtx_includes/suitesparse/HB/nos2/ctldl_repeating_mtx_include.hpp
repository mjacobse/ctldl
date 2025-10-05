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
  constexpr auto diag = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsityStatic<3, 0>();
  constexpr auto outer = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeStart{diag, next, outer, permutation};
}

constexpr auto getRepeatingMtxSparsityTridiag() {
  constexpr auto diag =
      ctldl::makeSparsityStatic<3, 3>({{0, 0}, {1, 1}, {2, 2}});
  constexpr auto subdiag =
      ctldl::makeSparsityStatic<3, 3>({{0, 0}, {1, 1}, {2, 1}, {1, 2}, {2, 2}});
  constexpr auto permutation = ctldl::PermutationStatic<3>{{0, 1, 2}};
  return ctldl::SparsityToFactorizeTridiagonal{diag, subdiag, permutation};
}

constexpr auto getRepeatingMtxSparsityLink() {
  constexpr auto prev = ctldl::makeEmptySparsityStatic<0, 3>();
  constexpr auto diag = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeLink{prev, diag, next, permutation};
}

constexpr auto getRepeatingMtxSparsityOuter() {
  constexpr auto subdiag = ctldl::makeEmptySparsityStatic<0, 3>();
  constexpr auto diag = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeOuter{subdiag, diag, permutation};
}
