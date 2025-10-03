#include <ctldl/factor_data/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/factorize/regularization_none.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <array>
#include <cstdio>
#include <vector>


namespace {

constexpr int dim = 10;

class MatrixA {
 public:
  static constexpr auto sparsity = ctldl::makeSparsity<dim, dim>({
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

  explicit MatrixA(const bool is_last_block) : m_is_last_block(is_last_block) {}

  constexpr double valueAt(const std::size_t i) const {
    if (m_is_last_block) {
      return values_last_block[i];
    }
    return values[i];
  }

 private:
  bool m_is_last_block;
};

struct MatrixB {
  static constexpr auto sparsity = ctldl::makeSparsity<dim, dim>({
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
      -0.2,
      -0.071554176, -0.035777088, -0.071554176, 0.035777088,
      -0.035777088, -0.017888544, -0.2, 0.035777088, -0.017888544,
      // empty row
      -0.2,
      -0.071554176, -0.035777088, -0.071554176, 0.035777088,
       -0.035777088, -0.017888544, -0.2, 0.035777088, -0.017888544,
      // empty row
      -0.2}};

  static constexpr double valueAt(const std::size_t i) {
    return values[i];
  }
};

constexpr ctldl::PermutationStatic<dim> permutation{{7, 8, 0, 4, 3, 2, 6, 5, 9, 1}};

}  // anonymous namespace


int main() {
  const int num_repetitions = 9;

  ctldl::FactorizationRepeatingBlockTridiagonal<
      MatrixA::sparsity, MatrixB::sparsity, double, permutation>
      factorization(num_repetitions);

  const auto matrix_values_A = [] {
    std::vector<MatrixA> ret(num_repetitions + 1, MatrixA(false));
    ret.back() = MatrixA(true);
    return ret;
  }();
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
