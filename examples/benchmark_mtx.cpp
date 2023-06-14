#include "ctldl_repeating_mtx_include.hpp"

#include <ctldl/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/fileio/mtx_file_read_repeating_block_tridiagonal.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <span>
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

  constexpr double valueAt(const std::size_t i) const { return values[i]; }
};

using MatrixA = MatrixInput<SparsityA>;
using MatrixB = MatrixInput<SparsityB>;
using Factorization =
    ctldl::FactorizationRepeatingBlockTridiagonal<SparsityA, SparsityB, double,
                                                  Permutation>;
using RightHandSide = std::vector<std::array<double, dim>>;


[[gnu::noinline]] void runFactorize(
    const std::span<const MatrixA> matrix_values_A,
    const std::span<const MatrixB> matrix_values_B,
    Factorization& factorization) {
  factorization.factor(matrix_values_A, matrix_values_B);
}

void benchmarkFactorize(benchmark::State& state,
                        const std::span<const MatrixA> matrix_values_A,
                        const std::span<const MatrixB> matrix_values_B) {
  const auto num_repetitions = std::size_t{matrix_values_B.size()};
  Factorization factorization(num_repetitions);

  auto matrix_values_A_data = matrix_values_A.data();
  auto matrix_values_B_data = matrix_values_B.data();
  auto factorization_A_data = factorization.dataA();
  auto factorization_B_data = factorization.dataB();
  for (auto _ : state) {
    benchmark::DoNotOptimize(matrix_values_A_data);
    benchmark::DoNotOptimize(matrix_values_B_data);
    benchmark::DoNotOptimize(factorization_A_data);
    benchmark::DoNotOptimize(factorization_B_data);
    benchmark::ClobberMemory();
    runFactorize(matrix_values_A, matrix_values_B, factorization);
    benchmark::ClobberMemory();
  }
}

[[gnu::noinline]] void runSolve(const Factorization& factorization,
                                RightHandSide& rhs_in_solution_out) {
  factorization.solveInPlace(rhs_in_solution_out);
}

void benchmarkSolve(benchmark::State& state,
                    const std::span<const MatrixA> matrix_values_A,
                    const std::span<const MatrixB> matrix_values_B) {
  const auto num_repetitions = std::size_t{matrix_values_B.size()};

  Factorization factorization(num_repetitions);
  factorization.factor(matrix_values_A, matrix_values_B);

  const auto rhs_single = [] {
    std::array<double, dim> rhs;
    std::fill(rhs.begin(), rhs.end(), 1.0);
    return rhs;
  }();
  const RightHandSide rhs(num_repetitions + 1, rhs_single);
  auto rhs_in_solution_out = rhs;

  auto factorization_A_data = factorization.dataA();
  auto factorization_B_data = factorization.dataB();
  auto solution_data = rhs_in_solution_out.data();
  for (auto _ : state) {
    benchmark::DoNotOptimize(factorization_A_data);
    benchmark::DoNotOptimize(factorization_B_data);
    benchmark::DoNotOptimize(solution_data);
    std::copy(rhs.cbegin(), rhs.cend(), rhs_in_solution_out.begin());
    benchmark::ClobberMemory();
    runSolve(factorization, rhs_in_solution_out);
    benchmark::ClobberMemory();
  }
}

[[gnu::noinline]] void runCombined(
    const std::span<const MatrixA> matrix_values_A,
    const std::span<const MatrixB> matrix_values_B,
    Factorization& factorization, RightHandSide& rhs_in_solution_out) {
  factorization.factor(matrix_values_A, matrix_values_B);
  factorization.solveInPlace(rhs_in_solution_out);
}

void benchmarkCombined(benchmark::State& state,
                       const std::span<const MatrixA> matrix_values_A,
                       const std::span<const MatrixB> matrix_values_B) {
  const auto num_repetitions = std::size_t{matrix_values_B.size()};
  Factorization factorization(num_repetitions);

  const auto rhs_single = [] {
    std::array<double, dim> rhs;
    std::fill(rhs.begin(), rhs.end(), 1.0);
    return rhs;
  }();
  const RightHandSide rhs(num_repetitions + 1, rhs_single);
  auto rhs_in_solution_out = rhs;

  auto matrix_values_A_data = matrix_values_A.data();
  auto matrix_values_B_data = matrix_values_B.data();
  auto factorization_A_data = factorization.dataA();
  auto factorization_B_data = factorization.dataB();
  auto solution_data = rhs_in_solution_out.data();
  for (auto _ : state) {
    benchmark::DoNotOptimize(matrix_values_A_data);
    benchmark::DoNotOptimize(matrix_values_B_data);
    benchmark::DoNotOptimize(factorization_A_data);
    benchmark::DoNotOptimize(factorization_B_data);
    benchmark::DoNotOptimize(solution_data);
    std::copy(rhs.cbegin(), rhs.cend(), rhs_in_solution_out.begin());
    benchmark::ClobberMemory();
    runCombined(matrix_values_A, matrix_values_B, factorization,
                rhs_in_solution_out);
    benchmark::ClobberMemory();
  }
}

}  // anonymous namespace


int main(int argc, char** argv) {
  if (argc < 2) {
    std::puts("Usage: benchmark_* FILEPATH_MTX\n");
    return EXIT_FAILURE;
  }

  const char* path_mtx = argv[1];
  const auto [matrix_values_A, matrix_values_B] =
      ctldl::mtxFileReadRepeatingBlockTridiagonal<MatrixA, MatrixB>(path_mtx);

  const auto unit = benchmark::kMicrosecond;
  benchmark::RegisterBenchmark("factorize", benchmarkFactorize, matrix_values_A,
                               matrix_values_B)
      ->Unit(unit);
  benchmark::RegisterBenchmark("solve", benchmarkSolve, matrix_values_A,
                               matrix_values_B)
      ->Unit(unit);
  benchmark::RegisterBenchmark("combined", benchmarkCombined, matrix_values_A,
                               matrix_values_B)
      ->Unit(unit);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
