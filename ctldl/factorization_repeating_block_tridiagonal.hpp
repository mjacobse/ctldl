#pragma once

#include <ctldl/factorization.hpp>
#include <ctldl/factorization_subdiagonal_block.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/solve_backward_substitution.hpp>
#include <ctldl/solve_forward_substitution.hpp>
#include <ctldl/sparsity/get_entries.hpp>
#include <ctldl/sparsity/is_nonzero_info.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/compute_elimination_tree_repeating.hpp>

#include <cstddef>
#include <memory>


namespace ctldl {

template <std::size_t dim>
struct IsNonZeroPair {
  IsNonzeroInfo<dim, dim> A;
  IsNonzeroInfo<dim, dim> B;
};

template <class Sparsity, class PermutationIn>
struct RepeatedSparsity {
  using SparsityInA = typename Sparsity::A;
  using SparsityInB = typename Sparsity::B;
  static_assert(SparsityInA::num_rows == SparsityInA::num_cols);
  static_assert(SparsityInB::num_cols == SparsityInA::num_cols);
  static_assert(SparsityInB::num_rows == SparsityInA::num_rows);
  static constexpr auto dim = std::size_t{SparsityInA::num_rows};
  static constexpr Permutation<dim> permutation{PermutationIn::permutation};

  using SparsityA =
      SparsityCSR<SparsityLowerTriangle<SparsityInA, PermutationIn>>;
  using SparsityB =
      SparsityCSR<SparsityPermuted<SparsityInB, PermutationIn, PermutationIn>>;
  static constexpr auto tree =
      computeEliminationTreeRepeating<SparsityA, SparsityB>();

  static constexpr auto is_nonzero_pair = [] {
    IsNonzeroInfo<dim, dim> is_nonzero_A;
    IsNonzeroInfo<dim, dim> is_nonzero_B;
    const auto mark_nonzero = [&is_nonzero_A, &is_nonzero_B](
                                  const std::size_t i, const std::size_t j) {
      if (j >= dim) {
        is_nonzero_A[i][j - dim] = true;
      } else {
        is_nonzero_B[i][j] = true;
      }
    };

    std::array<std::size_t, 2*dim> visitor;
    std::iota(visitor.begin(), visitor.end(), 0);
    for (std::size_t i = 0; i < dim; ++i) {
      constexpr auto row_offset = dim;
      constexpr auto col_offset = dim;
      foreachAncestorInSubtree<SparsityB>(tree, i, visitor, mark_nonzero,
                                          row_offset);
      foreachAncestorInSubtree<SparsityA>(tree, i, visitor, mark_nonzero,
                                          row_offset, col_offset);
    }
    return IsNonZeroPair<dim>{is_nonzero_A, is_nonzero_B};
  }();

  struct A {
    static constexpr auto num_rows = dim;
    static constexpr auto num_cols = dim;
    static constexpr auto entries =
        getEntries([] { return is_nonzero_pair.A; });
  };
  struct B {
    static constexpr auto num_rows = dim;
    static constexpr auto num_cols = dim;
    static constexpr auto entries =
        getEntries([] { return is_nonzero_pair.B; });
  };
};

template <class FactorData>
struct FactorSubdiagonalBlockData {
  using Sparsity = typename FactorData::Sparsity;
  using Value = typename FactorData::Value;

  template <class FactorDataAbove>
  FactorSubdiagonalBlockData(const FactorData& self,
                             const FactorDataAbove& above)
      : L(self.L), D(above.D) {
    static_assert(Sparsity::num_cols == FactorDataAbove::Sparsity::num_cols);
  }

  const std::array<Value, Sparsity::nnz>& L;
  const std::array<Value, Sparsity::num_cols>& D;
};

// Factorization of a matrix of the form
//
// [A  B'         ]
// [B  A  B'      ]
// [   :  :  :    ]
// [      B  A  B']
// [         B  A ]
template <class Sparsity, class Value,
          class PermutationIn = PermutationIdentity>
class FactorizationRepeatingBlockTridiagonal {
 private:
  using SparsityFactor = RepeatedSparsity<Sparsity, PermutationIn>;
  using SparsityFactorA = typename SparsityFactor::A;
  using SparsityFactorB = typename SparsityFactor::B;
  using FactorA = Factorization<SparsityFactorA, Value, PermutationIn>;
  using FactorB = FactorizationSubdiagonalBlock<SparsityFactorB, Value,
                                                PermutationIn, PermutationIn>;
  std::size_t m_num_repetitions;
  std::unique_ptr<FactorA[]> m_diag;
  std::unique_ptr<FactorB[]> m_subdiag;

 public:
  explicit FactorizationRepeatingBlockTridiagonal(
      const std::size_t num_repetitions)
      : m_num_repetitions(num_repetitions),
        m_diag(new FactorA[num_repetitions + 1]),
        m_subdiag(new FactorB[num_repetitions]) {}

  template <class MatrixValuesA, class MatrixValuesB>
  [[gnu::noinline]] void factor(const MatrixValuesA& values_A,
                                const MatrixValuesB& values_B) {
    m_diag[0].factor(values_A[0]);
    for (std::size_t i = 0; i < m_num_repetitions; ++i) {
      m_subdiag[i].factor(values_B[i], m_diag[i]);
      m_diag[i + 1].factor(values_A[i + 1],
                           FactorSubdiagonalBlockData{m_subdiag[i], m_diag[i]});
    }
  }

  template <class Rhs>
  [[gnu::noinline]] void solveInPlace(Rhs& rhs) const {
    forwardSolve(rhs);
    diagonalSolve(rhs);
    backwardSolve(rhs);
  }

  template <class Rhs>
  [[gnu::noinline]] void forwardSolve(Rhs& rhs) const {
    solveForwardSubstitution(m_diag[0], rhs[0]);
    for (std::size_t i = 1; i <= m_num_repetitions; ++i) {
      solveForwardSubstitution(m_diag[i], rhs[i], m_subdiag[i - 1], rhs[i - 1]);
    }
  }

  template <class Rhs>
  [[gnu::noinline]] void diagonalSolve(Rhs& rhs) const {
    for (std::size_t i = 0; i <= m_num_repetitions; ++i) {
      m_diag[i].diagonalSolve(rhs[i]);
    }
  }

  template <class Rhs>
  [[gnu::noinline]] void backwardSolve(Rhs& rhs) const {
    for (std::size_t i = m_num_repetitions; i > 0; --i) {
      solveBackwardSubstitution(m_diag[i], rhs[i], m_subdiag[i - 1], rhs[i - 1]);
    }
    solveBackwardSubstitution(m_diag[0], rhs[0]);
  }
};

}  // namespace ctldl
