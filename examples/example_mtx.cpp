#include <ctldl/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>

#include "ctldl_repeating_mtx_include.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
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

struct SparsityAtDiscretePoint {
  using A = SparsityA;
  using B = SparsityB;
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

using MatrixA = MatrixInput<SparsityAtDiscretePoint::A>;
using MatrixB = MatrixInput<SparsityAtDiscretePoint::B>;


class ReadMtxException : public std::runtime_error {
 public:
  explicit ReadMtxException(const char* error_message)
      : std::runtime_error(error_message) {}
};

void checkMtx(const bool condition, const char* message) {
  if (!condition) {
    throw ReadMtxException(message);
  }
}

std::pair<std::vector<MatrixA>, std::vector<MatrixB>> readRepeatingMtxFile(
    const char* filepath) {
  std::ifstream file(filepath);
  std::string line;

  while (std::getline(file, line) && line[0] == '%') {
    // skip comment lines
  }

  std::size_t num_repetitions;
  {
    std::size_t num_rows;
    std::size_t num_cols;
    std::size_t nnz;
    const auto num_read = std::sscanf(line.c_str(), " %zu %zu %zu", &num_rows,
                                      &num_cols, &nnz);
    checkMtx(num_read == 3, "Error reading matrix dimensions");
    checkMtx(num_rows == num_cols, "Matrix must be square");
    checkMtx(num_rows > 0, "Matrix dimension must be non-zero");
    checkMtx(num_rows % dim == 0,
             "Matrix dimension must be multiple of expected repeating block "
             "dimension");
    num_repetitions = (num_rows / dim) - 1;
  }

  std::vector<std::array<double, MatrixA::nnz>> values_A(num_repetitions + 1,
                                                         {0.0});
  std::vector<std::array<double, MatrixB::nnz>> values_B(num_repetitions,
                                                         {0.0});
  while (std::getline(file, line)) {
    std::size_t orig_row_index;
    std::size_t orig_col_index;
    double value;
    const auto num_read = std::sscanf(line.c_str(), " %zu %zu %lf",
                                      &orig_row_index, &orig_col_index, &value);
    checkMtx(num_read == 3, "Error reading matrix entry");
    checkMtx(orig_row_index >= orig_col_index,
             "Matrix entries must be in lower triangle");
    checkMtx(orig_row_index > 0, "Row index must be positive");
    checkMtx(orig_col_index > 0, "Col index must be negative");
    orig_row_index -= 1;
    orig_col_index -= 1;

    const std::size_t repetition_index = (orig_col_index / dim);
    const bool is_diagonal_block =
        (orig_row_index / dim == orig_col_index / dim);

    const auto row_index = orig_row_index % dim;
    const auto col_index = orig_col_index % dim;
    if (is_diagonal_block) {
      checkMtx(MatrixA::Sparsity::isNonZero(row_index, col_index),
               "Entry is not covered by compiled diagonal block sparsity");
      const auto entry_index =
          MatrixA::Sparsity::entryIndex(row_index, col_index);
      values_A[repetition_index][entry_index] = value;
    } else {
      checkMtx(MatrixB::Sparsity::isNonZero(row_index, col_index),
               "Entry is not covered by compiled subdiagonal block sparsity");
      const auto entry_index =
          MatrixB::Sparsity::entryIndex(row_index, col_index);
      values_B[repetition_index][entry_index] = value;
    }
  }

  std::vector<MatrixA> matrices_A(values_A.size(), MatrixA({0.0}));
  std::transform(values_A.cbegin(), values_A.cend(), matrices_A.begin(),
                 [](const auto& values) { return MatrixA{values}; });
  std::vector<MatrixB> matrices_B(values_B.size(), MatrixB({0.0}));
  std::transform(values_B.cbegin(), values_B.cend(), matrices_B.begin(),
                 [](const auto& values) { return MatrixB{values}; });

  return {matrices_A, matrices_B};
}

}  // anonymous namespace


int main(const int argc, const char** argv) {
  if (argc < 2 || argc > 3) {
    std::puts("Usage: example_* FILEPATH_MTX [NUM_ITERATIONS]\n");
    return EXIT_FAILURE;
  }

  const char* path_mtx = argv[1];
  const std::size_t num_iterations = (argc > 2) ? std::stoull(argv[2]) : 1;

  const auto [matrix_values_A, matrix_values_B] =
      readRepeatingMtxFile(path_mtx);
  const auto num_repetitions = matrix_values_B.size();

  ctldl::FactorizationRepeatingBlockTridiagonal<SparsityAtDiscretePoint, double,
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
    factorization.factor(matrix_values_A, matrix_values_B);
    std::copy(rhs.cbegin(), rhs.cend(), rhs_in_solution_out.begin());
    factorization.solveInPlace(rhs_in_solution_out);
  }

  for (const auto& solution_part : rhs_in_solution_out) {
    for (const auto v : solution_part) {
      std::printf("%f\n", v);
    }
  }
}
