#pragma once

#include <ctldl/factorization.hpp>
#include <ctldl/factorization_subdiagonal_block.hpp>
#include <ctldl/factorize/factorize_entry_wise.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/factorize/factorize_up_looking.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_identity.hpp>
#include <ctldl/solve_backward_substitution.hpp>
#include <ctldl/solve_forward_substitution.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/foreach_nonzero_with_fill_repeated.hpp>

#include <array>
#include <cstddef>
#include <memory>
#include <utility>


namespace ctldl {

template <class SparsityInA, class SparsityInB, class PermutationIn>
struct RepeatedSparsity {
  static_assert(SparsityInA::num_rows == SparsityInA::num_cols);
  static_assert(SparsityInB::num_cols == SparsityInA::num_cols);
  static_assert(SparsityInB::num_rows == SparsityInA::num_rows);
  static constexpr auto dim = std::size_t{SparsityInA::num_rows};
  static constexpr Permutation<dim> permutation{PermutationIn::permutation};

  using SparsityA =
      SparsityCSR<SparsityLowerTriangle<SparsityInA, PermutationIn>>;
  using SparsityB =
      SparsityCSR<SparsityPermuted<SparsityInB, PermutationIn, PermutationIn>>;

  static constexpr auto nnz_pair = [] {
    std::size_t nnz_A = 0;
    std::size_t nnz_B = 0;
    const auto count_nonzero = [&](const std::size_t /*i*/,
                                   const std::size_t j) {
      if (j >= dim) {
        nnz_A += 1;
      } else {
        nnz_B += 1;
      }
    };
    foreachNonZeroWithFillRepeated<SparsityA, SparsityB>(count_nonzero);
    return std::pair{nnz_A, nnz_B};
  }();

  static constexpr auto entries_pair = [] {
    std::array<Entry, nnz_pair.first> entries_A;
    std::array<Entry, nnz_pair.second> entries_B;
    std::size_t entry_index_A = 0;
    std::size_t entry_index_B = 0;
    const auto add_nonzero = [&](const std::size_t i, const std::size_t j) {
      if (j >= dim) {
        entries_A[entry_index_A] = Entry{i, j - dim};
        entry_index_A += 1;
      } else {
        entries_B[entry_index_B] = Entry{i, j};
        entry_index_B += 1;
      }
    };
    foreachNonZeroWithFillRepeated<SparsityA, SparsityB>(add_nonzero);
    // sorting is not needed for correctness, but helps performance
    sortEntriesRowMajorSorted(entries_A);
    sortEntriesRowMajorSorted(entries_B);
    return std::pair{entries_A, entries_B};
  }();

  struct A {
    static constexpr auto num_rows = dim;
    static constexpr auto num_cols = dim;
    static constexpr auto entries = entries_pair.first;
  };
  struct B {
    static constexpr auto num_rows = dim;
    static constexpr auto num_cols = dim;
    static constexpr auto entries = entries_pair.second;
  };
};

// Factorization of a matrix of the form
//
// [A  B'         ]
// [B  A  B'      ]
// [   :  :  :    ]
// [      B  A  B']
// [         B  A ]
template <class SparsityA, class SparsityB, class Value,
          class PermutationIn = PermutationIdentity>
class FactorizationRepeatingBlockTridiagonal {
 private:
  using SparsityFactor = RepeatedSparsity<SparsityA, SparsityB, PermutationIn>;
  using SparsityFactorA = typename SparsityFactor::A;
  using SparsityFactorB = typename SparsityFactor::B;
  using FactorA = Factorization<SparsityFactorA, Value, PermutationIn>;
  using FactorB = FactorizationSubdiagonalBlock<SparsityFactorB, Value,
                                                PermutationIn, PermutationIn>;
  std::size_t m_num_repetitions;
  std::unique_ptr<FactorA[]> m_diag;
  std::unique_ptr<FactorB[]> m_subdiag;

  template <class MatrixValuesA, class MatrixValuesB>
  void factorizeBlockRow(const FactorA& diag_previous,
                         const MatrixValuesB& values_B,
                         const MatrixValuesA& values_A, FactorB& subdiag,
                         FactorA& diag, FactorizeMethodEntryWise) {
    factorizeEntryWise(diag_previous, values_B, values_A, subdiag, diag);
  }

  template <class MatrixValuesA, class MatrixValuesB>
  void factorizeBlockRow(const FactorA& diag_previous,
                         const MatrixValuesB& values_B,
                         const MatrixValuesA& values_A, FactorB& subdiag,
                         FactorA& diag, FactorizeMethodUpLooking) {
    factorizeUpLooking(diag_previous, values_B, values_A, subdiag, diag);
  }

 public:
  explicit FactorizationRepeatingBlockTridiagonal(
      const std::size_t num_repetitions)
      : m_num_repetitions(num_repetitions),
        m_diag(new FactorA[num_repetitions + 1]),
        m_subdiag(new FactorB[num_repetitions]) {}

  const FactorA* dataA() const noexcept { return m_diag.get(); }
  const FactorB* dataB() const noexcept { return m_subdiag.get(); }

  template <class FactorizeMethodTag = FactorizeMethodUpLooking,
            class MatrixValuesA, class MatrixValuesB>
  void factor(const MatrixValuesA& values_A, const MatrixValuesB& values_B,
              const FactorizeMethodTag method_tag = {}) {
    m_diag[0].factor(values_A[0], method_tag);
    for (std::size_t i = 0; i < m_num_repetitions; ++i) {
      factorizeBlockRow(m_diag[i], values_B[i], values_A[i + 1], m_subdiag[i],
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
  void forwardSolve(Rhs& rhs) const {
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
  void backwardSolve(Rhs& rhs) const {
    solveBackwardSubstitution(m_diag[m_num_repetitions], rhs[m_num_repetitions]);
    for (std::size_t i = m_num_repetitions; i > 0; --i) {
      solveBackwardSubstitution(m_diag[i - 1], rhs[i - 1], m_subdiag[i - 1], rhs[i]);
    }
  }
};

}  // namespace ctldl
