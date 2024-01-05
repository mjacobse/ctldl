#pragma once

#include <ctldl/permutation/permutation_enumeration.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

template <class Value>
struct TestMatrixNos2A {
  struct Matrix {
    struct Entry {
      std::size_t row_index;
      std::size_t col_index;
      Value value;
    };
    struct Sparsity {
      static constexpr int num_rows = 3;
      static constexpr int num_cols = 3;
      static constexpr std::array<Entry, 3> entries = {
          {{0, 0, 640000.0}, {1, 1, 25600000.0}, {2, 2, 78643200000.0}}};
    };
    static constexpr auto sparsity = Sparsity{};

    static constexpr double valueAt(const std::size_t i) {
      return Sparsity::entries[i].value;
    }
  };

  static constexpr const char* description() { return "Nos2 A"; }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

template <class Value>
struct TestMatrixNos2B {
  struct Matrix {
    struct Entry {
      std::size_t row_index;
      std::size_t col_index;
      Value value;
    };
    struct Sparsity {
      static constexpr int num_rows = 3;
      static constexpr int num_cols = 3;
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

  static constexpr const char* description() { return "Nos2 B"; }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& /*value_generator*/,
                       const std::size_t num_matrices) {
    return std::vector<Matrix>(num_matrices);
  }
};

}  // namespace ctldl
