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
  constexpr auto diag = ctldl::makeSparsityStatic<4, 4>({
    {0,0},
    {1,1},
    {2,2},
    {3,3},
  });
  constexpr auto next = ctldl::makeSparsityStatic<7, 4>({
    {0,0},
    {0,1},
    {1,1},
    {2,2},
    {1,3},
    {2,3},
  });
  constexpr auto outer = ctldl::makeEmptySparsityStatic<0, 4>();
  constexpr auto permutation = ctldl::PermutationStatic<4>{{
    0,
    1,
    2,
    3,
  }};
  return ctldl::SparsityToFactorizeStart{diag, next, outer, permutation};
}

constexpr auto getRepeatingMtxSparsityTridiag() {
  constexpr auto diag = ctldl::makeSparsityStatic<7, 7>({
    {0,0},
    {3,0},
    {4,0},
    {1,1},
    {4,1},
    {6,1},
    {2,2},
    {5,2},
    {6,2},
    {3,3},
    {4,4},
    {5,5},
    {6,6},
  });
  constexpr auto subdiag = ctldl::makeSparsityStatic<7, 7>({
    {0,3},
    {0,4},
    {1,4},
    {2,5},
    {1,6},
    {2,6},
  });
  constexpr auto permutation = ctldl::PermutationStatic<7>{{
    3,
    5,
    4,
    0,
    1,
    2,
    6,
  }};
  return ctldl::SparsityToFactorizeTridiagonal{diag, subdiag, permutation};
}

constexpr auto getRepeatingMtxSparsityLink() {
  constexpr auto prev = ctldl::makeEmptySparsityStatic<0, 7>();
  constexpr auto diag = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeLink{prev, diag, next, permutation};
}

constexpr auto getRepeatingMtxSparsityOuter() {
  constexpr auto subdiag = ctldl::makeEmptySparsityStatic<0, 7>();
  constexpr auto diag = ctldl::makeEmptySparsityStatic<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeOuter{subdiag, diag, permutation};
}
