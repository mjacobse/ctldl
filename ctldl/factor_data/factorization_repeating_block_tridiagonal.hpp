#pragma once

#include <ctldl/factor_data/factorization_already_permuted.hpp>
#include <ctldl/factor_data/factorization_subdiagonal_block.hpp>
#include <ctldl/factorize/factorize_entry_wise.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/factorize/factorize_up_looking.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/solve/solve_backward_substitution.hpp>
#include <ctldl/solve/solve_forward_substitution.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/filled_in_sparsity_repeating.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <utility>


namespace ctldl {

// Factorization of a matrix of the form
//
// [A  B'         ]
// [B  A  B'      ]
// [   :  :  :    ]
// [      B  A  B']
// [         B  A ]
template <auto sparsity_in_A, auto sparsity_in_B, class Value_,
          auto permutation_in = PermutationIdentity{}>
class FactorizationRepeatingBlockTridiagonal {
 private:
  static constexpr auto sparsity_A = makeSparsity(sparsity_in_A);
  static constexpr auto sparsity_B = makeSparsity(sparsity_in_B);
  static_assert(isSquare(sparsity_A));
  static_assert(isSquare(sparsity_B));
  static_assert(sparsity_B.num_rows == sparsity_A.num_rows);

 public:
  using Value = Value_;
  static constexpr auto dim = std::size_t{sparsity_A.num_rows};
  static constexpr Permutation<dim> permutation{permutation_in};

  static constexpr auto sparsity_factor =
      getFilledInSparsityRepeating<sparsity_A, sparsity_B, permutation>();
  using FactorA = FactorizationAlreadyPermuted<sparsity_factor.diagonal, Value,
                                               permutation>;
  using FactorB =
      FactorizationSubdiagonalBlock<sparsity_factor.subdiagonal, Value,
                                    permutation, permutation>;

  explicit FactorizationRepeatingBlockTridiagonal(
      const std::size_t num_repetitions)
      : m_num_repetitions(num_repetitions),
        m_diag(new FactorA[num_repetitions + 1]),
        m_subdiag(new FactorB[num_repetitions]) {}

  std::span<const FactorA> blocksA() const noexcept {
    return {m_diag.get(), m_num_repetitions + 1};
  }
  std::span<const FactorB> blocksB() const noexcept {
    return {m_subdiag.get(), m_num_repetitions};
  }

  template <class FactorizeMethodTag = FactorizeMethodUpLooking,
            class MatrixValuesA, class MatrixValuesB>
  void factorize(const MatrixValuesA& values_A, const MatrixValuesB& values_B,
                 const FactorizeMethodTag method_tag = {}) {
    m_diag[0].factorize(values_A[0], method_tag);
    for (std::size_t i = 0; i < m_num_repetitions; ++i) {
      ::ctldl::factorize(m_diag[i], values_B[i], values_A[i + 1], m_subdiag[i],
                         m_diag[i + 1], method_tag);
    }
  }

  template <class Rhs>
  void solveInPlace(Rhs& rhs) const {
    forwardSolve(rhs);
    diagonalSolve(rhs);
    backwardSolve(rhs);
  }

  template <class Rhs>
  [[gnu::flatten]] void forwardSolve(Rhs& rhs) const {
    solveForwardSubstitution(m_diag[0], rhs[0]);
    for (std::size_t i = 1; i <= m_num_repetitions; ++i) {
      solveForwardSubstitution(m_diag[i], rhs[i], m_subdiag[i - 1], rhs[i - 1]);
    }
  }

  template <class Rhs>
  void diagonalSolve(Rhs& rhs) const {
    for (std::size_t i = 0; i <= m_num_repetitions; ++i) {
      m_diag[i].diagonalSolve(rhs[i]);
    }
  }

  template <class Rhs>
  [[gnu::flatten]] void backwardSolve(Rhs& rhs) const {
    solveBackwardSubstitution(m_diag[m_num_repetitions], rhs[m_num_repetitions]);
    for (std::size_t i = m_num_repetitions; i > 0; --i) {
      solveBackwardSubstitution(m_diag[i - 1], rhs[i - 1], m_subdiag[i - 1], rhs[i]);
    }
  }

 private:
  std::size_t m_num_repetitions;
  std::unique_ptr<FactorA[]> m_diag;
  std::unique_ptr<FactorB[]> m_subdiag;
};

}  // namespace ctldl
