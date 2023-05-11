#pragma once

#include <ctldl/factorize.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/filled_in_sparsity.hpp>
#include <ctldl/utility/square.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

template <class OriginalSparsity, class Value_,
          class PermutationIn = PermutationIdentity>
class Factorization {
 public:
  static_assert(OriginalSparsity::num_rows == OriginalSparsity::num_cols);
  static constexpr auto dim = std::size_t{OriginalSparsity::num_rows};
  using Sparsity = SparsityCSR<FilledInSparsity<OriginalSparsity>>;
  static constexpr auto nnz = std::size_t{Sparsity::nnz};
  using Value = Value_;
  static constexpr Permutation<dim> permutation{PermutationIn::permutation};
  static constexpr auto permutation_row = permutation;
  static constexpr auto permutation_col = permutation;

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
      const auto i_orig = permutation[i];
      rhs_in_solution_out[i_orig] /= D[i];
    }
  }
};

}  // namespace ctldl
