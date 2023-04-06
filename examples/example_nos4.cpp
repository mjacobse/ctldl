#include <ctldl/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstdio>
#include <vector>


namespace {

constexpr int dim = 10;

struct SparsityA {
  static constexpr int num_rows = dim;
  static constexpr int num_cols = dim;
  static constexpr int nnz = 16;
  static constexpr std::array<ctldl::Entry, nnz> entries = {{
    {0, 0},
    {1, 0}, {1, 1},
    {2, 0}, {2, 2},
    {3, 3},
    {4, 2}, {4, 4},
    {5, 5},
    {6, 4}, {6, 6},
    {7, 7},
    {8, 6}, {8, 8},
    {9, 8}, {9, 9}}};
};

struct SparsityB {
  static constexpr int num_rows = dim;
  static constexpr int num_cols = dim;
  static constexpr int nnz = 21;
  static constexpr std::array<ctldl::Entry, nnz> entries = {{
      {1, 1},
      {2, 0}, {2, 1}, {2, 4}, {2, 5},
      {3, 0}, {3, 1}, {3, 3}, {3, 4}, {3, 5},
      // empty row
      {5, 5},
      {6, 4}, {6, 5}, {6, 8}, {6, 9},
      {7, 4}, {7, 5}, {7, 7}, {7, 8}, {7, 9},
      // empty row
      {9, 9}}};
};

struct SparsityAtDiscretePoint {
  using A = SparsityA;
  using B = SparsityB;
};

template <class Sparsity_>
class Matrix {
 public:
  using Sparsity = ctldl::SparsityCSR<Sparsity_>;
  static constexpr std::size_t nnz = Sparsity::nnz;

  explicit Matrix(std::array<double, nnz> values)
      : m_values(values) {}

  constexpr double get(const std::size_t i, const std::size_t j) const {
    if (!Sparsity::is_nonzero[i][j]) {
      return 0.0;
    }
    return m_values[Sparsity::entryIndex(i, j)];
  }

 private:
  std::array<double, nnz> m_values;
};

struct Permutation {
  static constexpr std::array<std::size_t, dim> permutation = {7, 8, 0, 4, 3,
                                                               2, 6, 5, 9, 1};
};

}  // anonymous namespace


int main() {
  using MatrixA = Matrix<SparsityAtDiscretePoint::A>;
  using MatrixB = Matrix<SparsityAtDiscretePoint::B>;
  const int num_repetitions = 9;

  ctldl::FactorizationRepeatingBlockTridiagonal<SparsityAtDiscretePoint, double,
                                                Permutation>
      factorization(num_repetitions);

  const auto matrix_values_A = [&] {
    const MatrixA A({
      0.17155418,
      0.035777088, 0.41788854,
      -0.1, 0.34310835,
      0.43577709,
      -0.1, 0.34310835,
      0.43577709,
      -0.1, 0.34310835,
      0.43577709,
      -0.1, 0.17155418,
      -0.035777088, 0.41788854
    });
    const MatrixA A_last({
      0.1,
      0.0, 0.2,
      -0.1, 0.34310835,
      0.23577709,
      -0.1, 0.2,
      0.2,
      -0.1, 0.34310835,
      0.23577709,
      -0.1, 0.1,
      0.0, 0.2
    });
    std::vector<MatrixA> ret(num_repetitions + 1, A);
    ret.back() = A_last;
    return ret;
  }();
  const MatrixB B({
    -0.2,
    -0.071554176, -0.035777088, -0.071554176, 0.035777088,
    -0.035777088, -0.017888544, -0.2, 0.035777088, -0.017888544,
    // empty row
    -0.2,
    -0.071554176, -0.035777088, -0.071554176, 0.035777088,
     -0.035777088, -0.017888544, -0.2, 0.035777088, -0.017888544,
    // empty row
    -0.2
  });
  const std::vector<MatrixB> matrix_values_B(num_repetitions, B);

  for (int i = 0; i < 1000000; ++i) {
    factorization.factor(matrix_values_A, matrix_values_B);
  }

  const auto rhs_single = [] {
    std::array<double, dim> rhs;
    std::fill(rhs.begin(), rhs.end(), 1.0);
    return rhs;
  }();
  const std::vector<std::array<double, dim>> rhs(num_repetitions + 1,
                                                 rhs_single);
  auto rhs_in_solution_out = rhs;
  for (int i = 0; i < 1000000; ++i) {
    std::copy(rhs.cbegin(), rhs.cend(), rhs_in_solution_out.begin());
    factorization.solveInPlace(rhs_in_solution_out);
  }

  for (const auto& solution_part : rhs_in_solution_out) {
    for (const auto v : solution_part) {
      std::printf("%f\n", v);
    }
  }
}
