#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

template <class Value>
struct TestMatrixNos4A {
  static constexpr int block_dim = 10;

  class Matrix {
   public:
    static constexpr auto sparsity = makeSparsity<block_dim, block_dim>({
        {0, 0},
        {1, 0}, {1, 1},
        {2, 0}, {2, 2},
        {3, 3},
        {4, 2}, {4, 4},
        {5, 5},
        {6, 4}, {6, 6},
        {7, 7},
        {8, 6}, {8, 8},
        {9, 8}, {9, 9}});

    static constexpr std::array<double, sparsity.nnz()> values = {{
        0.17155418,
        0.035777088, 0.41788854,
        -0.1, 0.34310835,
        0.43577709,
        -0.1, 0.34310835,
        0.43577709,
        -0.1, 0.34310835,
        0.43577709,
        -0.1, 0.17155418,
        -0.035777088, 0.41788854}};
    static constexpr std::array<double, sparsity.nnz()> values_last_block = {{
        0.1,
        0.0, 0.2,
        -0.1, 0.34310835,
        0.23577709,
        -0.1, 0.2,
        0.2,
        -0.1, 0.34310835,
        0.23577709,
        -0.1, 0.1,
        0.0, 0.2}};

    explicit Matrix(const bool is_last_block)
        : m_is_last_block(is_last_block) {}

    constexpr auto valueAt(const std::size_t i) const {
      if (m_is_last_block) {
        return static_cast<Value>(values_last_block[i]);
      }
      return static_cast<Value>(values[i]);
    }

   private:
    bool m_is_last_block;
  };

  static constexpr const char* description() { return "Nos4 A"; }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    std::vector<Matrix> matrices(num_matrices, Matrix(false));
    if (!matrices.empty()) {
      matrices.back() = Matrix(true);
    }
    return matrices;
  }
};

template <class Value>
struct TestMatrixNos4B {
  static constexpr int block_dim = 10;
  struct Matrix {
    static constexpr auto sparsity = makeSparsity<block_dim, block_dim>({
        // empty row
        {1, 1},
        {2, 0}, {2, 1}, {2, 4}, {2, 5},
        {3, 0}, {3, 1}, {3, 3}, {3, 4}, {3, 5},
        // empty row
        {5, 5},
        {6, 4}, {6, 5}, {6, 8}, {6, 9},
        {7, 4}, {7, 5}, {7, 7}, {7, 8}, {7, 9},
        // empty row
        {9, 9}});

    static constexpr std::array<double, sparsity.nnz()> values = {{
        // empty row
        -0.2,
        -0.071554176, -0.035777088, -0.071554176, 0.035777088,
        -0.035777088, -0.017888544, -0.2, 0.035777088, -0.017888544,
        // empty row
        -0.2, -0.071554176, -0.035777088, -0.071554176, 0.035777088,
        -0.035777088, -0.017888544, -0.2, 0.035777088, -0.017888544,
        // empty row
        -0.2}};

    static constexpr auto valueAt(const std::size_t i) {
      return static_cast<Value>(values[i]);
    }
  };

  static constexpr const char* description() { return "Nos4 B"; }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

struct TestPermutationNos4 {
  static constexpr std::array<std::size_t, TestMatrixNos4A<double>::block_dim>
      value = {7, 8, 0, 4, 3, 2, 6, 5, 9, 1};
};

}  // namespace ctldl
