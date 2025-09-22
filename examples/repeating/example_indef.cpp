#include <ctldl/factor_data/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/factorize/regularization_small_positive_constant.hpp>
#include <ctldl/permutation/permutation.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <vector>


namespace {

struct Entry {
  std::size_t row_index;
  std::size_t col_index;
  double value;
};

constexpr int dim = 2;

struct MatrixA {
  struct Sparsity {
    static constexpr int num_rows = dim;
    static constexpr int num_cols = dim;
    static constexpr std::array<Entry, 1> entries = {{{1, 0, 1.0}}};
  };
  static constexpr auto sparsity = Sparsity{};

  static constexpr double valueAt(const std::size_t i) {
    return Sparsity::entries[i].value;
  }
};

struct MatrixB {
  struct Sparsity {
    static constexpr int num_rows = dim;
    static constexpr int num_cols = dim;
    static constexpr std::array<Entry, 5> entries = {};
  };
  static constexpr auto sparsity = Sparsity{};

  static constexpr double valueAt(const std::size_t i) {
    return Sparsity::entries[i].value;
  }
};

constexpr ctldl::Permutation<dim> permutation{{0, 1}};

}  // anonymous namespace


int main() {
  const int num_repetitions = 1;

  ctldl::FactorizationRepeatingBlockTridiagonal<
      MatrixA::sparsity, MatrixB::sparsity, float, permutation>
      factorization(num_repetitions);

  const std::vector<MatrixA> matrix_values_A(num_repetitions + 1);
  const std::vector<MatrixB> matrix_values_B(num_repetitions);

  factorization.factorize(matrix_values_A, matrix_values_B,
                          ctldl::RegularizationSmallPositiveConstant{});

  std::vector<std::array<float, dim>> rhs = {{1.0, 2.0}, {3.0, 4.0}};
  factorization.solveInPlace(rhs);

  for (const auto& solution_part : rhs) {
    for (const auto v : solution_part) {
      std::printf("%f\n", static_cast<double>(v));
    }
  }
}
