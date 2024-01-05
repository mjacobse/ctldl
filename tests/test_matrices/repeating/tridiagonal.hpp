#pragma once

#include <ctldl/permutation/permutation_enumeration.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

template <int block_dim_, class Value>
struct TestMatrixTridiagonal {
  static constexpr int block_dim = block_dim_;

  struct Matrix {
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
        std::array<Entry, std::size_t{nnz}> ret;
        fixInitIfZeroLengthArray(ret);
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
    static constexpr auto sparsity = Sparsity{};

    constexpr Value valueAt(const std::size_t i) const {
      return Sparsity::entries[i].value;
    }
  };

  static constexpr const char* description() { return "Tridiagonal"; }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

template <int block_dim, class Value>
struct TestMatrixSingleEntryTopRight {
  struct Matrix {
    static constexpr auto sparsity =
        makeSparsity<block_dim, block_dim>({Entry{0, block_dim - 1}});

    constexpr auto valueAt(const std::size_t /*i*/) const {
      return Value{-1.0};
    }
  };

  static constexpr const char* description() {
    return "Single entry top right";
  }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

}  // namespace ctldl
