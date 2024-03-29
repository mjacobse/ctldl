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
  constexpr auto diag = ctldl::makeSparsity<11, 11>({
    {0,0},
    {1,1},
    {2,2},
    {5,2},
    {7,2},
    {3,3},
    {4,4},
    {5,5},
    {6,6},
    {7,7},
    {8,8},
    {9,9},
    {10,10},
  });
  constexpr auto next = ctldl::makeSparsity<20, 11>({
    {0,0},
    {0,1},
    {1,1},
    {2,2},
    {3,2},
    {2,3},
    {3,3},
    {1,4},
    {3,4},
    {4,4},
    {5,5},
    {5,6},
    {6,6},
    {6,7},
    {7,7},
    {8,8},
    {4,9},
    {8,9},
    {7,10},
    {8,10},
  });
  constexpr auto outer = ctldl::makeSparsity<1, 11>({
    {0,1},
    {0,2},
    {0,3},
    {0,4},
    {0,6},
    {0,7},
    {0,9},
    {0,10},
  });
  constexpr auto permutation = ctldl::Permutation<11>{{
    6,
    0,
    7,
    3,
    4,
    2,
    5,
    9,
    8,
    10,
    1,
  }};
  return ctldl::SparsityToFactorizeStart{diag, next, outer, permutation};
}

constexpr auto getRepeatingMtxSparsityTridiag() {
  constexpr auto diag = ctldl::makeSparsity<20, 20>({
    {9,0},
    {10,0},
    {10,1},
    {13,1},
    {11,2},
    {12,2},
    {11,3},
    {12,3},
    {13,3},
    {14,3},
    {16,3},
    {13,4},
    {18,4},
    {14,5},
    {15,5},
    {15,6},
    {16,6},
    {16,7},
    {19,7},
    {17,8},
    {18,8},
    {19,8},
    {9,9},
    {10,10},
    {11,11},
    {14,11},
    {16,11},
    {12,12},
    {13,13},
    {14,14},
    {16,14},
    {15,15},
    {16,16},
    {17,17},
    {18,18},
    {19,19},
  });
  constexpr auto subdiag = ctldl::makeSparsity<20, 20>({
    {0,9},
    {0,10},
    {1,10},
    {2,11},
    {3,11},
    {2,12},
    {3,12},
    {1,13},
    {3,13},
    {4,13},
    {3,14},
    {5,14},
    {5,15},
    {6,15},
    {3,16},
    {6,16},
    {7,16},
    {8,17},
    {4,18},
    {8,18},
    {7,19},
    {8,19},
  });
  constexpr auto permutation = ctldl::Permutation<20>{{
    11,
    12,
    13,
    18,
    9,
    15,
    14,
    10,
    19,
    17,
    1,
    16,
    0,
    4,
    8,
    3,
    6,
    7,
    5,
    2,
  }};
  return ctldl::SparsityToFactorizeTridiagonal{diag, subdiag, permutation};
}

constexpr auto getRepeatingMtxSparsityLink() {
  constexpr auto prev = ctldl::makeEmptySparsity<0, 20>();
  constexpr auto diag = ctldl::makeEmptySparsity<0, 0>();
  constexpr auto next = ctldl::makeEmptySparsity<1, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeLink{prev, diag, next, permutation};
}

constexpr auto getRepeatingMtxSparsityOuter() {
  constexpr auto subdiag = ctldl::makeSparsity<1, 20>({
    {0,0},
    {0,1},
    {0,2},
    {0,3},
    {0,4},
    {0,5},
    {0,6},
    {0,7},
    {0,8},
    {0,10},
    {0,11},
    {0,12},
    {0,13},
    {0,14},
    {0,15},
    {0,16},
    {0,18},
    {0,19},
  });
  constexpr auto diag = ctldl::makeSparsity<1, 1>({
    {0,0},
  });
  constexpr auto permutation = ctldl::Permutation<1>{{
    0,
  }};
  return ctldl::SparsityToFactorizeOuter{subdiag, diag, permutation};
}
