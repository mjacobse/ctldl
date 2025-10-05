#pragma once

#include <ctldl/factor_data/factorization_repeating_block_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factorize/regularization.hpp>
#include <ctldl/matrix/matrix_link.hpp>
#include <ctldl/matrix/matrix_outer.hpp>
#include <ctldl/matrix/matrix_start.hpp>
#include <ctldl/matrix/matrix_tridiagonal.hpp>
#include <ctldl/matrix/matrix_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/vector/vector_tridiagonal_arrowhead_linked.hpp>

#include <array>
#include <cstddef>


namespace ctldl {

// Factorization of a matrix of the form
//
// [A  B'         ]
// [B  A  B'      ]
// [   :  :  :    ]
// [      B  A  B']
// [         B  A ]
template <SparsityStatic sparsity_A, SparsityStatic sparsity_B, class Value_,
          PermutationStatic<sparsity_A.numRows()> permutation_in =
              PermutationIdentity{}>
class FactorizationRepeatingBlockTridiagonal {
 private:
  static_assert(isSquare(sparsity_A));
  static_assert(isSquare(sparsity_B));
  static_assert(sparsity_B.numRows() == sparsity_A.numRows());

  static constexpr auto sparsity_to_factorize =
      SparsityToFactorizeTridiagonalArrowheadLinked{
          makeEmptySparsityToFactorizeStart<0, sparsity_A.numRows(), 0>(),
          SparsityToFactorizeTridiagonal{sparsity_A, sparsity_B,
                                         permutation_in},
          makeEmptySparsityToFactorizeLink<sparsity_A.numCols(), 0, 0>(),
          makeEmptySparsityToFactorizeOuter<sparsity_A.numCols(), 0>()};
  using Base = FactorizationRepeatingBlockTridiagonalArrowheadLinked<
      sparsity_to_factorize, Value_>;
  Base m_base;

 public:
  using Value = Value_;
  static constexpr auto dim = std::size_t{sparsity_A.numRows()};

  explicit FactorizationRepeatingBlockTridiagonal(
      const std::size_t num_repetitions)
      : m_base(num_repetitions) {}

  const auto& data() const noexcept { return m_base.data().tridiag; }

  template <class FactorizeMethodTag = FactorizeMethodUpLooking,
            class MatrixValuesA, class MatrixValuesB>
  void factorize(const MatrixValuesA& values_A, const MatrixValuesB& values_B,
                 const Regularization auto& regularization,
                 const FactorizeMethodTag method_tag = {}) {
    m_base.factorize(
        MatrixTridiagonalArrowheadLinked{
            makeEmptyMatrixStart<0, dim, 0>(),
            MatrixTridiagonal{values_A, values_B},
            makeEmptyMatrixLink<dim, 0, 0>(),
            makeEmptyMatrixRepeatingOuter<0, dim>(m_base.numRepetitions())},
        regularization, method_tag);
  }

  template <class Rhs>
  void solveInPlace(Rhs& rhs) const {
    VectorTridiagonalArrowheadLinked<std::array<Value, 0>, Rhs&,
                                     std::array<Value, 0>, std::array<Value, 0>>
        rhs_extended{{}, rhs, {}, {}};
    m_base.solveInPlace(rhs_extended);
  }
};

}  // namespace ctldl
