#pragma once

#include <ctldl/factorize.hpp>
#include <ctldl/sparsity/filled_in_sparsity.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/utility/square.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

template <class OriginalSparsity, class Value_>
class Factorization {
 public:
  static_assert(OriginalSparsity::num_rows == OriginalSparsity::num_cols);
  static constexpr auto dim = std::size_t{OriginalSparsity::num_rows};
  using Sparsity = SparsityCSR<FilledInSparsity<OriginalSparsity>>;
  static constexpr auto nnz = std::size_t{Sparsity::nnz};
  using Value = Value_;

  template <class Matrix>
  void factor(const Matrix& matrix) {
    factorize(*this, matrix);
  }

  template <class Matrix, class FactorDataLeft>
  void factor(const Matrix& matrix, const FactorDataLeft& left) {
    factorize(*this, matrix, left);
  }

  template <class Rhs>
  void diagonalSolve(Rhs& rhs_in_solution_out) const {
    diagonalSolveImpl(rhs_in_solution_out.data());
  }

  std::array<Value, nnz> L;
  std::array<Value, dim> D;

 private:
  void diagonalSolveImpl(Value* __restrict rhs_in_solution_out) const {
    for (std::size_t i = 0; i < dim; ++i) {
      rhs_in_solution_out[i] /= D[i];
    }
  }
};

}  // namespace ctldl
