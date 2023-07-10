#pragma once

#include <ctldl/factorize/factorize.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/filled_in_sparsity.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <auto sparsity_orig_in, class Value_,
          auto permutation_in = PermutationIdentity{}>
class Factorization {
 public:
  static constexpr auto sparsity_orig = makeSparsity(sparsity_orig_in);
  static_assert(sparsity_orig.num_rows == sparsity_orig.num_cols);
  static constexpr auto dim = std::size_t{sparsity_orig.num_rows};
  using Value = Value_;
  static constexpr Permutation<dim> permutation{permutation_in};

  static constexpr auto sparsity =
      SparsityCSR(getFilledInSparsity<sparsity_orig>());
  static constexpr auto nnz = std::size_t{sparsity.nnz};
  static constexpr auto permutation_row = permutation;
  static constexpr auto permutation_col = permutation;

  template <class FactorizeMethodTag = FactorizeMethodUpLooking, class Matrix>
  void factorize(const Matrix& matrix,
                 const FactorizeMethodTag method_tag = {}) {
    ::ctldl::factorize(*this, matrix, method_tag);
  }

  template <class Rhs>
  void diagonalSolve(Rhs& rhs_in_solution_out) const {
    diagonalSolveImpl(rhs_in_solution_out.data());
  }

  std::array<Value, nnz> L;
  std::array<Value, dim> D;

 private:
  template <class ValueRhs>
  void diagonalSolveImpl(ValueRhs* __restrict rhs_in_solution_out) const {
    for (std::size_t i = 0; i < dim; ++i) {
      const auto i_orig = permutation[i];
      rhs_in_solution_out[i_orig] /= D[i];
    }
  }
};

}  // namespace ctldl
