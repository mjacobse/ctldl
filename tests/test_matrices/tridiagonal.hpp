#pragma once

#include <ctldl/permutation/permutation_enumeration.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

template <int block_dim_, class Value>
struct TestMatrixTridiagonal {
  static constexpr int block_dim = block_dim_;

  struct MatrixA {
    struct Entry {
      std::size_t row_index;
      std::size_t col_index;
      Value value;
    };

    struct Sparsity {
      static constexpr int num_rows = block_dim;
      static constexpr int num_cols = block_dim;
      static constexpr int nnz = block_dim + (block_dim - 1);
      static constexpr auto entries = [] {
        std::array<Entry, std::size_t{nnz}> ret{};
        std::size_t entry_index = 0;
        for (std::size_t i = 0; i < num_rows; ++i) {
          ret[entry_index] = {i, i, 2.0};
          entry_index += 1;
        }
        for (std::size_t i = 1; i < num_rows; ++i) {
          ret[entry_index] = {i, i - 1, -1.0};
          entry_index += 1;
        }
        return ret;
      }();
    };

    constexpr Value valueAt(const std::size_t i) const {
      return Sparsity::entries[i].value;
    }
  };

  struct MatrixB {
    struct Sparsity {
      static constexpr int num_rows = block_dim;
      static constexpr int num_cols = block_dim;
      static constexpr int nnz = 1;
      static constexpr std::array<ctldl::Entry, nnz> entries = {
          {{0, block_dim - 1}}};
    };

    constexpr auto valueAt(const std::size_t /*i*/) const {
      return Value{-1.0};
    }
  };

  static constexpr double expected_error_amplifier = 4096.0;
  static constexpr const char* description() { return "Tridiagonal"; }

  explicit TestMatrixTridiagonal(const std::size_t num_repetitions)
      : matrices_A(num_repetitions + 1), matrices_B(num_repetitions) {}

  std::vector<MatrixA> matrices_A;
  std::vector<MatrixB> matrices_B;
};

}  // namespace ctldl
