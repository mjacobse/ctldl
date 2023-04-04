#include <ctldl/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/sparsity/entry.hpp>

#include <array>
#include <cstdio>
#include <vector>


namespace {

constexpr int dim_A = 3;

template <int n>
struct SparsityDiagonalAndSubdiagonal {
  static constexpr int num_rows = n;
  static constexpr int num_cols = n;
  static constexpr int nnz = n + (n - 1);
  static constexpr auto entries = [] {
    std::array<ctldl::Entry, nnz> ret{};
    int entry_index = 0;
    for (int i = 0; i < num_rows; ++i) {
      ret[entry_index] = ctldl::Entry{i, i};
      entry_index += 1;
    }
    for (int i = 1; i < num_rows; ++i) {
      ret[entry_index] = ctldl::Entry{i, i - 1};
      entry_index += 1;
    }
    return ret;
  }();
};

struct SparsityB {
  static constexpr int num_rows = dim_A;
  static constexpr int num_cols = dim_A;
  static constexpr int nnz = 1;
  static constexpr std::array<ctldl::Entry, nnz> entries = {{{0, dim_A - 1}}};
};

struct SparsityAtDiscretePoint {
  using A = SparsityDiagonalAndSubdiagonal<dim_A>;
  using B = SparsityB;
};

template <class Sparsity_>
class Matrix {
 public:
  using Sparsity = ctldl::SparsityCSR<Sparsity_>;
  static constexpr int nnz = Sparsity::nnz;

  explicit Matrix(std::array<double, nnz> values)
      : m_values(values) {}

  constexpr double get(int i, int j) const {
    if (!Sparsity::is_nonzero[i][j]) {
      return 0.0;
    }
    return m_values[Sparsity::entryIndex(i, j)];
  }

 private:
  std::array<double, nnz> m_values;
};

}  // anonymous namespace


int main() {
  using MatrixA = Matrix<SparsityAtDiscretePoint::A>;
  using MatrixB = Matrix<SparsityAtDiscretePoint::B>;
  const int num_repetitions = 20 - 1;

  ctldl::FactorizationRepeatingBlockTridiagonal<SparsityAtDiscretePoint, double>
      factorization(num_repetitions);

  const auto values_A = [] {
    std::array<double, MatrixA::nnz> values;
    std::fill_n(values.begin(), dim_A, 2.0);
    std::fill(values.begin() + dim_A, values.end(), -1.0);
    values[0] = 2.0;
    for (int i = 1; i < MatrixA::nnz; i += 2) {
      values[i] = -1.0;
      values[i + 1] = 2.0;
    }
    return values;
  }();
  const MatrixA A{values_A};
  const MatrixB B{{-1.0}};
  const std::vector<MatrixA> matrix_values_A(num_repetitions + 1, A);
  const std::vector<MatrixB> matrix_values_B(num_repetitions, B);

  for (int i = 0; i < 10000000; ++i) {
    factorization.factor(matrix_values_A, matrix_values_B);
  }

  const auto rhs_single = [] {
    std::array<double, dim_A> rhs;
    std::fill(rhs.begin(), rhs.end(), 1.0);
    return rhs;
  }();
  const std::vector<std::array<double, dim_A>> rhs(num_repetitions + 1,
                                                   rhs_single);
  auto rhs_in_solution_out = rhs;
  for (int i = 0; i < 10000000; ++i) {
    std::copy(rhs.cbegin(), rhs.cend(), rhs_in_solution_out.begin());
    factorization.solveInPlace(rhs_in_solution_out);
  }

  for (const auto solution_part : rhs_in_solution_out) {
    for (const auto v : solution_part) {
      std::printf("%f\n", v);
    }
  }
}
