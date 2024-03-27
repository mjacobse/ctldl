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
  constexpr auto next = ctldl::makeEmptySparsity<22, 0>();
  constexpr auto outer = ctldl::makeEmptySparsity<2, 0>();
  constexpr ctldl::Permutation<0> permutation{};
  return ctldl::SparsityToFactorizeStart{diag, next, outer, permutation};
}

constexpr auto getRepeatingMtxSparsityTridiag() {
  constexpr auto diag = ctldl::makeSparsity<22, 22>({
    {0,0},
    {1,0},
    {3,0},
    {4,0},
    {8,0},
    {9,0},
    {10,0},
    {11,0},
    {12,0},
    {13,0},
    {14,0},
    {1,1},
    {3,1},
    {4,1},
    {8,1},
    {9,1},
    {10,1},
    {11,1},
    {12,1},
    {13,1},
    {15,1},
    {2,2},
    {3,2},
    {7,2},
    {9,2},
    {14,2},
    {15,2},
    {16,2},
    {17,2},
    {3,3},
    {9,3},
    {10,3},
    {11,3},
    {12,3},
    {13,3},
    {14,3},
    {15,3},
    {17,3},
    {4,4},
    {7,4},
    {8,4},
    {9,4},
    {18,4},
    {21,4},
    {5,5},
    {6,5},
    {9,5},
    {19,5},
    {20,5},
    {6,6},
    {9,6},
    {16,6},
    {20,6},
    {7,7},
    {9,7},
    {17,7},
    {21,7},
    {8,8},
    {9,8},
    {18,8},
    {9,9},
    {19,9},
    {10,10},
    {11,11},
    {12,12},
    {13,13},
  });
  constexpr auto subdiag = ctldl::makeSparsity<22, 22>({
    {0,14},
    {2,14},
    {3,14},
    {1,15},
    {2,15},
    {3,15},
    {2,16},
    {6,16},
    {2,17},
    {3,17},
    {7,17},
    {4,18},
    {8,18},
    {5,19},
    {9,19},
    {5,20},
    {6,20},
    {4,21},
    {7,21},
  });
  constexpr auto permutation = ctldl::Permutation<22>{{
    18,
    10,
    17,
    21,
    13,
    12,
    11,
    20,
    16,
    8,
    14,
    15,
    4,
    7,
    19,
    1,
    3,
    2,
    0,
    6,
    9,
    5,
  }};
  return ctldl::SparsityToFactorizeTridiagonal{diag, subdiag, permutation};
}

constexpr auto getRepeatingMtxSparsityLink() {
  constexpr auto prev = ctldl::makeSparsity<14, 22>({
    {0,14},
    {2,14},
    {3,14},
    {1,15},
    {2,15},
    {3,15},
    {2,16},
    {6,16},
    {2,17},
    {3,17},
    {7,17},
    {4,18},
    {8,18},
    {5,19},
    {9,19},
    {5,20},
    {6,20},
    {4,21},
    {7,21},
  });
  constexpr auto diag = ctldl::makeSparsity<14, 14>({
    {0,0},
    {10,0},
    {11,0},
    {12,0},
    {13,0},
    {1,1},
    {10,1},
    {11,1},
    {12,1},
    {13,1},
    {2,2},
    {3,2},
    {7,2},
    {3,3},
    {10,3},
    {11,3},
    {12,3},
    {13,3},
    {4,4},
    {7,4},
    {5,5},
    {6,5},
    {6,6},
    {7,7},
    {8,8},
    {9,9},
    {10,10},
    {11,11},
    {12,12},
    {13,13},
  });
  constexpr auto next = ctldl::makeSparsity<2, 14>({
    {0,2},
    {0,3},
    {0,4},
    {0,5},
    {0,6},
    {0,7},
    {0,8},
    {0,9},
  });
  constexpr auto permutation = ctldl::Permutation<14>{{
    2,
    4,
    10,
    9,
    5,
    11,
    13,
    8,
    12,
    3,
    7,
    6,
    1,
    0,
  }};
  return ctldl::SparsityToFactorizeLink{prev, diag, next, permutation};
}

constexpr auto getRepeatingMtxSparsityOuter() {
  constexpr auto subdiag = ctldl::makeSparsity<2, 22>({
    {0,0},
    {0,1},
    {0,2},
    {0,3},
    {0,4},
    {0,5},
    {0,6},
    {0,7},
    {0,8},
    {0,9},
    {1,13},
    {0,14},
    {1,14},
    {0,15},
    {1,15},
    {0,16},
    {1,16},
    {0,17},
    {1,17},
    {0,18},
    {1,18},
    {0,19},
    {1,19},
    {0,20},
    {1,20},
    {0,21},
    {1,21},
  });
  constexpr auto diag = ctldl::makeSparsity<2, 2>({
    {0,0},
    {1,1},
  });
  constexpr auto permutation = ctldl::Permutation<2>{{
    1,
    0,
  }};
  return ctldl::SparsityToFactorizeOuter{subdiag, diag, permutation};
}
