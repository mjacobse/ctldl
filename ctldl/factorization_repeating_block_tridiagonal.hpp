#pragma once

#include <ctldl/factorization.hpp>
#include <ctldl/factorization_subdiagonal_block.hpp>
#include <ctldl/sparsity/get_entries.hpp>
#include <ctldl/sparsity/get_filled_in_is_nonzero_info.hpp>
#include <ctldl/sparsity/get_is_nonzero_info.hpp>
#include <ctldl/sparsity/is_nonzero_info.hpp>

#include <memory>


namespace ctldl {

template <int dim>
struct IsNonZeroPair {
  IsNonzeroInfo<dim, dim> A;
  IsNonzeroInfo<dim, dim> B;
};

template <class Sparsity>
struct RepeatedSparsity {
  static constexpr auto is_nonzero_pair = [] {
    static_assert(Sparsity::A::num_rows == Sparsity::A::num_cols);
    static_assert(Sparsity::B::num_cols == Sparsity::A::num_cols);
    static_assert(Sparsity::B::num_rows == Sparsity::A::num_rows);
    constexpr int dim = Sparsity::A::num_rows;

    auto is_nonzero_A = getIsNonzeroInfo<typename Sparsity::A>();
    auto is_nonzero_B = getIsNonzeroInfo<typename Sparsity::B>();

    for (int iter = 0; iter < dim + 1; ++iter) {
      is_nonzero_A = getFilledInIsNonzeroInfo(is_nonzero_A);

      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          for (int k = 0; k < j; ++k) {
            if (is_nonzero_B[i][k] && is_nonzero_A[j][k]) {
              is_nonzero_B[i][j] = true;
            }
          }
        }
      }

      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < i; ++j) {
          for (int k = 0; k < dim; ++k) {
            if (is_nonzero_B[i][k] && is_nonzero_B[j][k]) {
              is_nonzero_A[i][j] = true;
            }
          }
        }
      }
    }
    return IsNonZeroPair<dim>{is_nonzero_A, is_nonzero_B};
  }();

  struct A {
    static constexpr int num_rows = Sparsity::A::num_rows;
    static constexpr int num_cols = Sparsity::A::num_cols;
    static constexpr auto entries =
        getEntries([] { return is_nonzero_pair.A; });
  };
  struct B {
    static constexpr int num_rows = Sparsity::B::num_rows;
    static constexpr int num_cols = Sparsity::B::num_cols;
    static constexpr auto entries =
        getEntries([] { return is_nonzero_pair.B; });
  };
};

// Factorization of a matrix of the form
//
// [A  B'         ]
// [B  A  B'      ]
// [   :  :  :    ]
// [      B  A  B']
// [         B  A ]
template <class Sparsity, class Value>
class FactorizationRepeatingBlockTridiagonal {
 private:
  using SparsityFactor = RepeatedSparsity<Sparsity>;
  using SparsityFactorA = typename SparsityFactor::A;
  using SparsityFactorB = typename SparsityFactor::B;
  using FactorA = Factorization<SparsityFactorA, Value>;
  using FactorB = FactorizationSubdiagonalBlock<SparsityFactorB, Value>;
  int m_num_repetitions;
  std::unique_ptr<FactorA[]> m_diag;
  std::unique_ptr<FactorB[]> m_subdiag;

 public:
  explicit FactorizationRepeatingBlockTridiagonal(const int num_repetitions)
      : m_num_repetitions(num_repetitions),
        m_diag(new FactorA[num_repetitions + 1]),
        m_subdiag(new FactorB[num_repetitions]) {
  }

  template <class MatrixValuesA, class MatrixValuesB>
  [[gnu::noinline]] void factor(const MatrixValuesA& values_A,
                                const MatrixValuesB& values_B) {
    m_diag[0].init(values_A[0]);
    m_diag[0].factor();
    for (int i = 0; i < m_num_repetitions; ++i) {
      m_subdiag[i].init(values_B[i]);
      m_subdiag[i].factor(m_diag[i]);
      m_diag[i + 1].init(values_A[i + 1]);
      m_subdiag[i].contribute(m_diag[i], m_diag[i + 1]);
      m_diag[i + 1].factor();
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
    m_diag[0].forwardSolve(rhs[0]);
    for (int i = 1; i <= m_num_repetitions; ++i) {
      m_subdiag[i - 1].forwardSolve(rhs[i - 1], rhs[i]);
      m_diag[i].forwardSolve(rhs[i]);
    }
  }

  template <class Rhs>
  [[gnu::noinline]] void diagonalSolve(Rhs& rhs) const {
    for (int i = 0; i <= m_num_repetitions; ++i) {
      m_diag[i].diagonalSolve(rhs[i]);
    }
  }

  template <class Rhs>
  [[gnu::noinline]] void backwardSolve(Rhs& rhs) const {
    for (int i = m_num_repetitions; i > 0; --i) {
      m_diag[i].backwardSolve(rhs[i]);
      m_subdiag[i - 1].backwardSolve(rhs[i - 1], rhs[i]);
    }
    m_diag[0].backwardSolve(rhs[0]);
  }
};

}  // namespace ctldl
