#pragma once

#include <ctldl/sparsity/sparsity.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class Value>
struct TestMatrixLFAT5 {
  static constexpr int dim = 14;
  struct Matrix {
    static constexpr auto sparsity = makeSparsity<dim, dim>(
        {{0, 0},   {3, 0},   {4, 0},   {1, 1},   {5, 1},   {2, 2},
         {6, 2},   {3, 3},   {7, 3},   {8, 3},   {4, 4},   {7, 4},
         {8, 4},   {5, 5},   {9, 5},   {6, 6},   {10, 6},  {7, 7},
         {11, 7},  {12, 7},  {8, 8},   {11, 8},  {12, 8},  {9, 9},
         {10, 10}, {11, 11}, {13, 11}, {12, 12}, {13, 12}, {13, 13}});

    static constexpr std::array<double, sparsity.nnz> values = {
        {1.57088000000000005e+00,  -9.42527999999999935e+01,
         7.85440000000000027e-01,  1.25664000000000000e+07,
         -6.28320000000000000e+06, 6.08806201550387560e-01,
         -3.04403100775193780e-01, 1.50804479999999967e+04,
         -7.54022399999999834e+03, 9.42527999999999935e+01,
         3.14176000000000011e+00,  -9.42527999999999935e+01,
         7.85440000000000027e-01,  1.25664000000000000e+07,
         -6.28320000000000000e+06, 6.08806201550387560e-01,
         -3.04403100775193780e-01, 1.50804479999999967e+04,
         -7.54022399999999834e+03, 9.42527999999999935e+01,
         3.14176000000000011e+00,  -9.42527999999999935e+01,
         7.85440000000000027e-01,  1.25664000000000000e+07,
         6.08806201550387560e-01,  1.50804479999999967e+04,
         9.42527999999999935e+01,  3.14176000000000011e+00,
         7.85440000000000027e-01,  1.57088000000000005e+00}};

    static constexpr auto valueAt(const std::size_t i) {
      return static_cast<Value>(values[i]);
    }
  };

  static constexpr double expected_error_amplifier = 8192.0;
  static constexpr const char* description() { return "LFAT5"; }

  static auto generate(const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

struct TestPermutationLFAT5 {
  static constexpr std::array<std::size_t, TestMatrixLFAT5<double>::dim>
      permutation = {13, 11, 12, 8, 7, 0, 4, 3, 9, 1, 5, 10, 2, 6};
};

}  // namespace ctldl
