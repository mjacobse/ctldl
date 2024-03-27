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
  constexpr auto next = ctldl::makeEmptySparsity<16, 0>();
  constexpr auto outer = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeStart{diag, next, outer, permutation};
}

constexpr auto getRepeatingMtxSparsityTridiag() {
  constexpr auto diag = ctldl::makeSparsity<16, 16>(
      {{0, 0},   {1, 0},   {1, 1},   {2, 2},   {3, 2},   {3, 3},   {4, 4},
       {5, 4},   {5, 5},   {6, 4},   {6, 5},   {6, 6},   {7, 4},   {7, 5},
       {7, 6},   {7, 7},   {8, 8},   {9, 8},   {9, 9},   {10, 10}, {11, 10},
       {11, 11}, {12, 12}, {13, 12}, {13, 13}, {14, 14}, {15, 14}, {15, 15}});
  constexpr auto subdiag = ctldl::makeSparsity<16, 16>(
      {{0, 1},   {1, 1},   {2, 2},   {3, 2},   {2, 3},   {3, 3},   {4, 5},
       {5, 5},   {6, 5},   {7, 5},   {4, 7},   {5, 7},   {6, 7},   {7, 7},
       {8, 9},   {9, 9},   {10, 11}, {11, 11}, {12, 12}, {13, 12}, {12, 13},
       {13, 13}, {14, 14}, {15, 14}, {14, 15}, {15, 15}});
  constexpr auto permutation = ctldl::Permutation<16>{
      {0, 1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15}};
  return ctldl::SparsityToFactorizeTridiagonal{diag, subdiag, permutation};
}

constexpr auto getRepeatingMtxSparsityLink() {
  constexpr auto prev = ctldl::makeEmptySparsity<0, 16>();
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeLink{prev, diag, next, permutation};
}

constexpr auto getRepeatingMtxSparsityOuter() {
  constexpr auto subdiag = ctldl::makeEmptySparsity<0, 16>();
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeOuter{subdiag, diag, permutation};
}
