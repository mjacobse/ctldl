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
  constexpr auto next = ctldl::makeEmptySparsity<10, 0>();
  constexpr auto outer = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeStart{diag, next, outer, permutation};
}

constexpr auto getRepeatingMtxSparsityTridiag() {
  constexpr auto diag = ctldl::makeSparsity<10, 10>({
    {0,0},
    {1,0},
    {2,0},
    {1,1},
    {2,2},
    {4,2},
    {3,3},
    {4,4},
    {6,4},
    {5,5},
    {6,6},
    {8,6},
    {7,7},
    {8,8},
    {9,8},
    {9,9},
  });
  constexpr auto subdiag = ctldl::makeSparsity<10, 10>({
    {2,0},
    {3,0},
    {1,1},
    {2,1},
    {3,1},
    {3,3},
    {2,4},
    {3,4},
    {6,4},
    {7,4},
    {2,5},
    {3,5},
    {5,5},
    {6,5},
    {7,5},
    {7,7},
    {6,8},
    {7,8},
    {6,9},
    {7,9},
    {9,9},
  });
  constexpr auto permutation =
      ctldl::PermutationStatic<10>{{7, 8, 0, 4, 3, 2, 6, 5, 9, 1}};
  return ctldl::SparsityToFactorizeTridiagonal{diag, subdiag, permutation};
}

constexpr auto getRepeatingMtxSparsityLink() {
  constexpr auto prev = ctldl::makeEmptySparsity<0, 10>();
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeLink{prev, diag, next, permutation};
}

constexpr auto getRepeatingMtxSparsityOuter() {
  constexpr auto subdiag = ctldl::makeEmptySparsity<0, 10>();
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::PermutationStatic<0> permutation{};
  return ctldl::SparsityToFactorizeOuter{subdiag, diag, permutation};
}
