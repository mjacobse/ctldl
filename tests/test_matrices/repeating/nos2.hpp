#pragma once

#include <ctldl/permutation/permutation_enumeration.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

template <class Value>
struct TestMatrixNos2 {
  struct Entry {
    std::size_t row_index;
    std::size_t col_index;
    Value value;
  };

  static constexpr int block_dim = 3;

  struct MatrixA {
    struct Sparsity {
      static constexpr int num_rows = block_dim;
      static constexpr int num_cols = block_dim;
      static constexpr std::array<Entry, 3> entries = {
          {{0, 0, 640000.0}, {1, 1, 25600000.0}, {2, 2, 78643200000.0}}};
    };
    static constexpr auto sparsity = Sparsity{};

    static constexpr double valueAt(const std::size_t i) {
      return Sparsity::entries[i].value;
    }
  };

  struct MatrixB {
    struct Sparsity {
      static constexpr int num_rows = block_dim;
      static constexpr int num_cols = block_dim;
      static constexpr std::array<Entry, 5> entries = {{{0, 0, -320000.0},
                                                        {1, 1, -39321600000.0},
                                                        {1, 2, -614400000.0},
                                                        {2, 1, 614400000.0},
                                                        {2, 2, 6400000.0}}};
    };
    static constexpr auto sparsity = Sparsity{};

    static constexpr Value valueAt(const std::size_t i) {
      return Sparsity::entries[i].value;
    }
  };

  static constexpr double expected_error_amplifier = 32768.0;
  static constexpr const char* description() { return "Nos2"; }

  explicit TestMatrixNos2(const std::size_t num_repetitions)
      : matrices_A(num_repetitions + 1), matrices_B(num_repetitions) {}

  std::vector<MatrixA> matrices_A;
  std::vector<MatrixB> matrices_B;
};

}  // namespace ctldl
