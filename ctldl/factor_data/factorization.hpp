#pragma once

#include <ctldl/factorize/factorize.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/factorize/regularization.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/permutation/permuted_entry_lower_triangle.hpp>
#include <ctldl/solve/solve_backward_substitution.hpp>
#include <ctldl/solve/solve_forward_substitution.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/symbolic/filled_in_sparsity.hpp>

#include <array>
#include <cstddef>
#include <type_traits>

namespace ctldl {

template <Sparsity sparsity_orig, class Value_,
          Permutation<sparsity_orig.num_rows> permutation_in =
              PermutationIdentity{}>
class Factorization {
 public:
  static_assert(isSquare(sparsity_orig));
  static constexpr auto dim = std::size_t{sparsity_orig.num_rows};
  using Value = Value_;
  static constexpr auto permutation = permutation_in;

  /**
   * Given an entry location in the factor, returns the location of the
   * corresponding entry in the original matrix, i.e. the location before the
   * permutation for the factorization was applied.
   *
   * With this special case of a symmetric factorization, the returned entry for
   * the original matrix will always be in the lower triangle.
   */
  static constexpr Entry origEntry(const Entry factor_entry) {
    return permutedEntryLowerTriangle(factor_entry, permutation);
  };
  /**
   * For a given row index in the factor, returns the corresponding row index in
   * the original matrix, i.e. the row index before the permutation for the
   * factorization was applied.
   */
  static constexpr auto origRowIndex(const std::size_t factor_row_index) {
    return std::size_t{permutation[factor_row_index]};
  }
  /**
   * For a given column index in the factor, returns the corresponding column
   * index in the original matrix, i.e. the column index before the permutation
   * for the factorization was applied.
   */
  static constexpr auto origColIndex(const std::size_t factor_col_index) {
    return std::size_t{permutation[factor_col_index]};
  }

  static constexpr auto sparsity =
      SparsityCSR(getFilledInSparsity<sparsity_orig, permutation>());
  static constexpr auto nnz = std::size_t{sparsity.nnz};
  static constexpr auto permutation_row = permutation;
  static constexpr auto permutation_col = permutation;

  template <class FactorizeMethodTag = FactorizeMethodUpLooking, class Matrix>
  void factorize(const Matrix& matrix,
                 const Regularization auto& regularization,
                 const FactorizeMethodTag method_tag = {}) {
    ::ctldl::factorize(*this, matrix, regularization, method_tag);
  }

  template <class Rhs>
  void solveInPlace(Rhs& rhs) const {
    solveForwardSubstitution(*this, rhs);
    diagonalSolve(rhs);
    solveBackwardSubstitution(*this, rhs);
  }

  template <class Rhs>
  void diagonalSolve(Rhs&& rhs_in_solution_out) const {
    using ValueRhs = std::remove_reference_t<decltype(rhs_in_solution_out[0])>;
    for (std::size_t i = 0; i < dim; ++i) {
      const auto i_orig = origRowIndex(i);
      rhs_in_solution_out[i_orig] /= static_cast<ValueRhs>(D[i]);
    }
  }

  std::array<Value, nnz> L;
  std::array<Value, dim> D;
};

}  // namespace ctldl
