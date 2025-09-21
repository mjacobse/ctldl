#include <ctldl/factor_data/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/factorize/regularization_none.hpp>
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

constexpr int dim = 3;

struct MatrixA {
  struct Sparsity {
    static constexpr int num_rows = dim;
    static constexpr int num_cols = dim;
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
    static constexpr int num_rows = dim;
    static constexpr int num_cols = dim;
    static constexpr std::array<Entry, 5> entries = {{{0, 0, -320000.0},
                                                      {1, 1, -39321600000.0},
                                                      {1, 2, -614400000.0},
                                                      {2, 1, 614400000.0},
                                                      {2, 2, 6400000.0}}};
  };
  static constexpr auto sparsity = Sparsity{};

  static constexpr double valueAt(const std::size_t i) {
    return Sparsity::entries[i].value;
  }
};

constexpr ctldl::Permutation<dim> permutation{{0, 1, 2}};

}  // anonymous namespace


int main() {
  const int num_repetitions = 318;

  ctldl::FactorizationRepeatingBlockTridiagonal<
      MatrixA::sparsity, MatrixB::sparsity, double, permutation>
      factorization(num_repetitions);

  const std::vector<MatrixA> matrix_values_A(num_repetitions + 1);
  const std::vector<MatrixB> matrix_values_B(num_repetitions);

  for (int i = 0; i < 1000000; ++i) {
    factorization.factorize(matrix_values_A, matrix_values_B,
                            ctldl::RegularizationNone{});
  }

  const auto rhs_single = [] {
    std::array<double, dim> rhs;
    std::ranges::fill(rhs, 1.0);
    return rhs;
  }();
  const std::vector<std::array<double, dim>> rhs(num_repetitions + 1,
                                                 rhs_single);
  auto rhs_in_solution_out = rhs;
  for (int i = 0; i < 1000000; ++i) {
    std::ranges::copy(rhs, rhs_in_solution_out.begin());
    factorization.solveInPlace(rhs_in_solution_out);
  }

  for (const auto& solution_part : rhs_in_solution_out) {
    for (const auto v : solution_part) {
      std::printf("%f\n", v);
    }
  }
}
