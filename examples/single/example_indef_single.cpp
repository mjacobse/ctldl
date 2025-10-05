#include <ctldl/factor_data/factorization.hpp>
#include <ctldl/factorize/regularization_small_positive_constant.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>


namespace {

constexpr int dim = 2;
constexpr ctldl::PermutationStatic<dim> permutation{{0, 1}};

struct Matrix {
  static constexpr auto sparsity = ctldl::makeSparsityStatic<dim, dim>({{1, 0}});

  std::array<double, sparsity.nnz()> values;
  constexpr double valueAt(const std::size_t i) const {
    return values[i];
  }
};

}  // anonymous namespace


int main() {
  ctldl::Factorization<Matrix::sparsity, float, permutation> factorization;

  Matrix matrix({1.0});
  factorization.factorize(matrix, ctldl::RegularizationSmallPositiveConstant{});
  std::array<float, dim> rhs = {1.0, 1.0};
  factorization.solveInPlace(rhs);

  for (const auto v : rhs) {
    std::printf("%f\n", static_cast<double>(v));
  }
}
