#include "ctldl_repeating_mtx_include.hpp"

#include <ctldl/factor_data/factorization_repeating_block_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factorize/regularization_small_positive_constant.hpp>
#include <ctldl/fileio/mtx_file_read_repeating_block_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/vector/vector_tridiagonal_arrowhead_linked.hpp>

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <vector>


namespace {

constexpr auto sparsity = ctldl::SparsityToFactorizeTridiagonalArrowheadLinked{
    getRepeatingMtxSparsityStart(),
    getRepeatingMtxSparsityTridiag(),
    getRepeatingMtxSparsityLink(),
    getRepeatingMtxSparsityOuter(),
};

using Factorization =
    ctldl::FactorizationRepeatingBlockTridiagonalArrowheadLinked<sparsity,
                                                                 double>;
using RightHandSide = ctldl::VectorTridiagonalArrowheadLinked<
    std::array<double, sparsity.dim_start>,
    std::vector<std::array<double, sparsity.dim_tridiag>>,
    std::array<double, sparsity.dim_link>,
    std::array<double, sparsity.dim_outer>>;
using Matrix =
    decltype(ctldl::mtxFileReadRepeatingBlockTridiagonalArrowheadLinked<
             sparsity>(""));

const auto factorize_method = ctldl::FactorizeMethodUpLooking{};

template <std::size_t dim>
auto getArrayOfOnes() {
  std::array<double, dim> arr;
  std::ranges::fill(arr, 1.0);
  return arr;
}

void copyRightHandSide(const RightHandSide& in, RightHandSide& out) {
  out.start = in.start;
  std::ranges::copy(in.tridiag, out.tridiag.begin());
  out.link = in.link;
  out.outer = in.outer;
}

[[gnu::noinline]] void runFactorize(const Matrix& matrix,
                                    Factorization& factorization) {
  factorization.factorize(matrix, ctldl::RegularizationSmallPositiveConstant{},
                          factorize_method);
}

void benchmarkFactorize(benchmark::State& state,
                        const Matrix& matrix) {
  const auto num_repetitions = std::size_t{matrix.tridiag.subdiag.size()};
  Factorization factorization(num_repetitions);

  auto matrix_values_A_data = matrix.tridiag.diag.data();
  auto matrix_values_B_data = matrix.tridiag.subdiag.data();
  auto factorization_A_data = factorization.data().tridiag.diag.data();
  auto factorization_B_data = factorization.data().tridiag.subdiag.data();
  for (auto _ : state) {
    benchmark::DoNotOptimize(matrix_values_A_data);
    benchmark::DoNotOptimize(matrix_values_B_data);
    benchmark::DoNotOptimize(factorization_A_data);
    benchmark::DoNotOptimize(factorization_B_data);
    benchmark::ClobberMemory();
    runFactorize(matrix, factorization);
    benchmark::ClobberMemory();
  }
}

[[gnu::noinline]] void runSolve(const Factorization& factorization,
                                RightHandSide& rhs_in_solution_out) {
  factorization.solveInPlace(rhs_in_solution_out);
}

void benchmarkSolve(benchmark::State& state, const Matrix& matrix) {
  const auto num_repetitions = std::size_t{matrix.tridiag.subdiag.size()};

  Factorization factorization(num_repetitions);
  factorization.factorize(matrix, ctldl::RegularizationSmallPositiveConstant{},
                          factorize_method);

  const RightHandSide rhs{
      getArrayOfOnes<sparsity.dim_start>(),
      {num_repetitions + 1, getArrayOfOnes<sparsity.dim_tridiag>()},
      getArrayOfOnes<sparsity.dim_link>(),
      getArrayOfOnes<sparsity.dim_outer>()};
  auto rhs_in_solution_out = rhs;

  auto factorization_A_data = factorization.data().tridiag.diag.data();
  auto factorization_B_data = factorization.data().tridiag.diag.data();
  auto solution_data = rhs_in_solution_out.tridiag.data();
  for (auto _ : state) {
    benchmark::DoNotOptimize(factorization_A_data);
    benchmark::DoNotOptimize(factorization_B_data);
    benchmark::DoNotOptimize(solution_data);
    copyRightHandSide(rhs, rhs_in_solution_out);
    benchmark::ClobberMemory();
    runSolve(factorization, rhs_in_solution_out);
    benchmark::ClobberMemory();
  }
}

[[gnu::noinline]] void runCombined(const Matrix& matrix,
                                   Factorization& factorization,
                                   RightHandSide& rhs_in_solution_out) {
  factorization.factorize(matrix, ctldl::RegularizationSmallPositiveConstant{},
                          factorize_method);
  factorization.solveInPlace(rhs_in_solution_out);
}

void benchmarkCombined(benchmark::State& state, const Matrix& matrix) {
  const auto num_repetitions = std::size_t{matrix.tridiag.subdiag.size()};
  Factorization factorization(num_repetitions);

  const RightHandSide rhs{
      getArrayOfOnes<sparsity.dim_start>(),
      {num_repetitions + 1, getArrayOfOnes<sparsity.dim_tridiag>()},
      getArrayOfOnes<sparsity.dim_link>(),
      getArrayOfOnes<sparsity.dim_outer>()};
  auto rhs_in_solution_out = rhs;

  auto matrix_values_A_data = matrix.tridiag.diag.data();
  auto matrix_values_B_data = matrix.tridiag.subdiag.data();
  auto factorization_A_data = factorization.data().tridiag.diag.data();
  auto factorization_B_data = factorization.data().tridiag.diag.data();
  auto solution_data = rhs_in_solution_out.tridiag.data();
  for (auto _ : state) {
    benchmark::DoNotOptimize(matrix_values_A_data);
    benchmark::DoNotOptimize(matrix_values_B_data);
    benchmark::DoNotOptimize(factorization_A_data);
    benchmark::DoNotOptimize(factorization_B_data);
    benchmark::DoNotOptimize(solution_data);
    copyRightHandSide(rhs, rhs_in_solution_out);
    benchmark::ClobberMemory();
    runCombined(matrix, factorization, rhs_in_solution_out);
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
  const auto matrix =
      ctldl::mtxFileReadRepeatingBlockTridiagonalArrowheadLinked<sparsity>(
          path_mtx);

  const auto unit = benchmark::kMicrosecond;
  benchmark::RegisterBenchmark("factorize", benchmarkFactorize, matrix)
      ->Unit(unit);
  benchmark::RegisterBenchmark("solve", benchmarkSolve, matrix)
      ->Unit(unit);
  benchmark::RegisterBenchmark("combined", benchmarkCombined, matrix)
      ->Unit(unit);

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
}
