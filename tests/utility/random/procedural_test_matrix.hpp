#pragma once

#include "tests/test_matrices/single/random.hpp"
#include "tests/utility/random/bool_matrix_distribution.hpp"
#include "tests/utility/random/compile_time_generate.hpp"
#include "tests/utility/random/linear_congruence_generator.hpp"

#include <ctldl/sparsity/bool_matrix.hpp>

#include <cstddef>

namespace ctldl {

enum class RandomMatrixSparsityKind {
  Any, Symmetric
};

template <RandomMatrixSparsityKind kind, std::size_t num_rows,
          std::size_t num_cols>
constexpr BoolMatrix<num_rows, num_cols> applyRandomMatrixSparsityKind(
    BoolMatrix<num_rows, num_cols> is_nonzero) {
  if constexpr (kind == RandomMatrixSparsityKind::Symmetric) {
    static_assert(num_rows == num_cols);
    for (std::size_t i = 0; i < num_rows; ++i) {
      for (std::size_t j = i + 1; j < num_cols; ++j) {
        is_nonzero.values[i][j] = false;
      }
    }
    for (std::size_t i = 0; i < num_rows; ++i) {
      is_nonzero.values[i][i] = true;
    }
  }
  return is_nonzero;
}

template <std::size_t num_rows, std::size_t num_cols,
          RandomMatrixSparsityKind kind, LinearCongruenceGenerator generator>
struct TestMatrixFromGenerator {
 private:
  using Value = double;
  static constexpr auto is_nonzero_generated = compileTimeGenerate(
      BoolMatrixDistribution<num_rows, num_cols>{}, generator);
  static constexpr auto is_nonzero =
      applyRandomMatrixSparsityKind<kind>(is_nonzero_generated.result);
  static constexpr auto sparsity = makeSparsity<num_rows, num_cols>(
      makeSparseEntriesFromBoolMatrix<is_nonzero>());
  using Base = TestMatrixRandom<sparsity, Value>;

 public:
  static constexpr auto next_generator = is_nonzero_generated.generator;
  using Matrix = typename Base::Matrix;

  static constexpr const char* description() { return Base::description(); }

  template <class ValueGenerator>
  static auto generate(ValueGenerator& value_generator,
                       const std::size_t num_matrices) {
    return Base::generate(value_generator, num_matrices);
  }
};

template <std::size_t seed, std::size_t num_rows, std::size_t num_cols,
          RandomMatrixSparsityKind kind = RandomMatrixSparsityKind::Any>
struct ProceduralTestMatrix {
  using type = TestMatrixFromGenerator<num_rows, num_cols, kind,
                                       LinearCongruenceGenerator{seed}>;
};

template <std::size_t seed, std::size_t num_rows, std::size_t num_cols>
using ProceduralTestMatrixT =
    typename ProceduralTestMatrix<seed, num_rows, num_cols>::type;

template <std::size_t seed, std::size_t dim>
using ProceduralTestMatrixSymmetricT =
    typename ProceduralTestMatrix<seed, dim, dim,
                                  RandomMatrixSparsityKind::Symmetric>::type;

}  // namespace ctldl
