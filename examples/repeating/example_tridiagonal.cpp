#include <ctldl/factor_data/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/factorize/regularization_none.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <array>
#include <cstdio>
#include <vector>


namespace {

constexpr int dim = 3;

struct MatrixA {
  struct Entry {
    std::size_t row_index;
    std::size_t col_index;
    double value;
  };

  struct Sparsity {
    static constexpr int num_rows = dim;
    static constexpr int num_cols = dim;
    static constexpr int nnz = dim + (dim - 1);
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
  static constexpr auto sparsity = Sparsity{};

  constexpr double valueAt(const std::size_t i) const {
    return Sparsity::entries[i].value;
  }
};

struct MatrixB {
  static constexpr auto sparsity =
      ctldl::makeSparsity<dim, dim>({ctldl::Entry{0, dim - 1}});

  constexpr double valueAt(const std::size_t /*i*/) const { return -1.0; }
};

}  // anonymous namespace


int main() {
  const int num_repetitions = 20 - 1;

  ctldl::FactorizationRepeatingBlockTridiagonal<MatrixA::sparsity,
                                                MatrixB::sparsity, double>
      factorization(num_repetitions);

  const std::vector<MatrixA> matrix_values_A(num_repetitions + 1);
  const std::vector<MatrixB> matrix_values_B(num_repetitions);
  for (int i = 0; i < 10000000; ++i) {
    factorization.factorize(matrix_values_A, matrix_values_B,
                            ctldl::RegularizationNone{});
  }

  const auto rhs_single = [] {
    std::array<double, dim> rhs;
    std::fill(rhs.begin(), rhs.end(), 1.0);
    return rhs;
  }();
  const std::vector<std::array<double, dim>> rhs(num_repetitions + 1,
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
