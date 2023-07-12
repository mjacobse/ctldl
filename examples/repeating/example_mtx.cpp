#include "ctldl_repeating_mtx_include.hpp"

#include <ctldl/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/fileio/mtx_file_read_repeating_block_tridiagonal.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <vector>


namespace {

constexpr int dim = getRepeatingMtxDim();
constexpr auto sparsity_A =
    ctldl::makeSparsity<dim, dim>(getRepeatingMtxEntriesA());
constexpr auto sparsity_B =
    ctldl::makeSparsity<dim, dim>(getRepeatingMtxEntriesB());
constexpr auto permutation = getRepeatingMtxPermutation();

template <auto sparsity_in>
struct MatrixInput {
  static constexpr auto sparsity = ctldl::SparsityCSR{sparsity_in};
  static constexpr auto nnz = std::size_t{sparsity.nnz};
  std::array<double, nnz> values;

  constexpr double valueAt(const std::size_t i) const {
    return values[i];
  }
};

using MatrixA = MatrixInput<sparsity_A>;
using MatrixB = MatrixInput<sparsity_B>;

}  // anonymous namespace


int main(const int argc, const char** argv) {
  if (argc < 2 || argc > 3) {
    std::puts("Usage: example_* FILEPATH_MTX [NUM_ITERATIONS]\n");
    return EXIT_FAILURE;
  }

  const char* path_mtx = argv[1];
  const std::size_t num_iterations = (argc > 2) ? std::stoull(argv[2]) : 1;

  const auto [matrix_values_A, matrix_values_B] =
      ctldl::mtxFileReadRepeatingBlockTridiagonal<MatrixA, MatrixB>(path_mtx);
  const auto num_repetitions = matrix_values_B.size();

  ctldl::FactorizationRepeatingBlockTridiagonal<sparsity_A, sparsity_B, double,
                                                permutation>
      factorization(num_repetitions);
  const auto rhs_single = [] {
    std::array<double, dim> rhs;
    std::fill(rhs.begin(), rhs.end(), 1.0);
    return rhs;
  }();
  const std::vector<std::array<double, dim>> rhs(num_repetitions + 1,
                                                 rhs_single);
  auto rhs_in_solution_out = rhs;

  for (std::size_t i = 0; i < num_iterations; ++i) {
    factorization.factorize(matrix_values_A, matrix_values_B);
    std::copy(rhs.cbegin(), rhs.cend(), rhs_in_solution_out.begin());
    factorization.solveInPlace(rhs_in_solution_out);
  }

  for (const auto& solution_part : rhs_in_solution_out) {
    for (const auto v : solution_part) {
      std::printf("%f\n", v);
    }
  }
}
