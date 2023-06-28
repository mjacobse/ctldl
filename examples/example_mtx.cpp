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

struct SparsityA {
  static constexpr int num_rows = dim;
  static constexpr int num_cols = dim;
  static constexpr auto entries = getRepeatingMtxEntriesA();
};

struct SparsityB {
  static constexpr int num_rows = dim;
  static constexpr int num_cols = dim;
  static constexpr auto entries = getRepeatingMtxEntriesB();
};

struct Permutation {
  static constexpr auto permutation = getRepeatingMtxPermutation();
};

template <class Sparsity_>
struct MatrixInput {
  using Sparsity = ctldl::SparsityCSR<Sparsity_>;
  static constexpr auto nnz = std::size_t{Sparsity::nnz};
  std::array<double, nnz> values;

  constexpr double valueAt(const std::size_t i) const {
    return values[i];
  }
};

using MatrixA = MatrixInput<SparsityA>;
using MatrixB = MatrixInput<SparsityB>;

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

  ctldl::FactorizationRepeatingBlockTridiagonal<SparsityA, SparsityB, double,
                                                Permutation>
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
