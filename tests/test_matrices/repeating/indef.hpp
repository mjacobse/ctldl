#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

template <class Value>
struct TestMatrixIndefA {
  struct Matrix {
    struct Entry {
      std::size_t row_index;
      std::size_t col_index;
      Value value;
    };
    struct Sparsity {
      static constexpr int num_rows = 2;
      static constexpr int num_cols = 2;
      static constexpr std::array<Entry, 1> entries = {{{1, 0, 1.0}}};
    };
    static constexpr auto sparsity = Sparsity{};

    static constexpr double valueAt(const std::size_t i) {
      return Sparsity::entries[i].value;
    }
  };

  static constexpr const char* description() { return "Indef A"; }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

template <class Value>
struct TestMatrixIndefB {
  struct Matrix {
    struct Entry {
      std::size_t row_index;
      std::size_t col_index;
      Value value;
    };
    struct Sparsity {
      static constexpr int num_rows = 2;
      static constexpr int num_cols = 2;
      static constexpr std::array<Entry, 0> entries = {};
    };
    static constexpr auto sparsity = Sparsity{};

    static constexpr Value valueAt(const std::size_t i) {
      return Sparsity::entries[i].value;
    }
  };

  static constexpr const char* description() { return "Indef B"; }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

}  // namespace ctldl
