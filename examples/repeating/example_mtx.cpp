#include "ctldl_repeating_mtx_include.hpp"

#include <ctldl/factor_data/factorization_repeating_block_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factorize/regularization_small_positive_constant.hpp>
#include <ctldl/fileio/mtx_file_read_repeating_block_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/vector/vector_tridiagonal_arrowhead_linked.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>


namespace {

constexpr auto sparsity = ctldl::SparsityToFactorizeTridiagonalArrowheadLinked{
    getRepeatingMtxSparsityStart(),
    getRepeatingMtxSparsityTridiag(),
    getRepeatingMtxSparsityLink(),
    getRepeatingMtxSparsityOuter(),
};

template <std::size_t dim>
auto getArrayOfOnes() {
  std::array<double, dim> arr;
  std::ranges::fill(arr, 1.0);
  return arr;
}

}  // anonymous namespace

int main(const int argc, const char** argv) {
  if (argc < 2 || argc > 3) {
    std::puts("Usage: example_* FILEPATH_MTX [NUM_ITERATIONS]\n");
    return EXIT_FAILURE;
  }

  const char* path_mtx = argv[1];
  const std::size_t num_iterations = (argc > 2) ? std::stoull(argv[2]) : 1;

  const auto matrix =
      ctldl::mtxFileReadRepeatingBlockTridiagonalArrowheadLinked<sparsity>(
          path_mtx);
  const auto num_repetitions = matrix.tridiag.subdiag.size();

  using Factorization =
      ctldl::FactorizationRepeatingBlockTridiagonalArrowheadLinked<sparsity,
                                                                   double>;
  Factorization factorization(num_repetitions);

  const ctldl::VectorTridiagonalArrowheadLinked rhs{
      getArrayOfOnes<sparsity.dim_start>(),
      std::vector{num_repetitions + 1, getArrayOfOnes<sparsity.dim_tridiag>()},
      getArrayOfOnes<sparsity.dim_link>(),
      getArrayOfOnes<sparsity.dim_outer>()};
  auto rhs_in_solution_out = rhs;

  for (std::size_t i = 0; i < num_iterations; ++i) {
    factorization.factorize(matrix,
                            ctldl::RegularizationSmallPositiveConstant{});

    rhs_in_solution_out.start = rhs.start;
    std::ranges::copy(rhs.tridiag, rhs_in_solution_out.tridiag.begin());
    rhs_in_solution_out.link = rhs.link;
    rhs_in_solution_out.outer = rhs.outer;

    factorization.solveInPlace(rhs_in_solution_out);
  }

  if (sparsity.dim_start > 0) {
    std::printf("------ start ------\n");
    for (const auto v : rhs_in_solution_out.start) {
      std::printf("%f\n", v);
    }
  }
  std::printf("----- tridiag -----\n");
  for (const auto& solution_part : rhs_in_solution_out.tridiag) {
    for (const auto v : solution_part) {
      std::printf("%f\n", v);
    }
  }
  if (sparsity.dim_link > 0) {
    std::printf("------ link -------\n");
    for (const auto v : rhs_in_solution_out.link) {
      std::printf("%f\n", v);
    }
  }
  if (sparsity.dim_outer > 0) {
    std::printf("------ outer ------\n");
    for (const auto v : rhs_in_solution_out.outer) {
      std::printf("%f\n", v);
    }
  }
}
