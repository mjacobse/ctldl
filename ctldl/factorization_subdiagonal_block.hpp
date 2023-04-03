#pragma once

#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/utility/square.hpp>

#include <array>

namespace ctldl {

template <class OriginalSparsity, class Value_>
class FactorizationSubdiagonalBlock {
 public:
  using Sparsity = SparsityCSR<OriginalSparsity>;
  static constexpr int nnz = Sparsity::nnz;
  static constexpr int num_rows = Sparsity::num_rows;
  using Value = Value_;

  template <class Matrix>
  void init(const Matrix& matrix) {
    // static_assert(isSparsitySubset(Matrix::Sparsity{}, Sparsity{}));
    initImplL(matrix);
  }

  template <class DiagonalBlock>
  void factor(const DiagonalBlock& above) {
    factorImpl(above);
  }

  template <class DiagonalBlock>
  void contribute(const DiagonalBlock& above, DiagonalBlock& right) const {
    contributeImpl(above, right);
    contributeImplDiag(above, right);
  }

  template <class Rhs>
  void forwardSolve(const Rhs& rhs, Rhs& solution) const {
    forwardSolveImpl(rhs.data(), solution.data());
  }

  template <class Rhs>
  void backwardSolve(Rhs& solution, const Rhs& rhs) const {
    backwardSolveImpl(solution.data(), rhs.data());
  }

 private:
  std::array<Value, Sparsity::nnz> L;

  template <int entry_index_ij = 0, class Matrix>
  [[gnu::always_inline]] void initImplL(const Matrix& matrix) {
    if constexpr (entry_index_ij < nnz) {
      constexpr auto i = Sparsity::entries[entry_index_ij].row_index;
      constexpr auto j = Sparsity::entries[entry_index_ij].col_index;
      if constexpr (Matrix::Sparsity::is_nonzero[i][j]) {
        L[entry_index_ij] = matrix.m_values[Matrix::Sparsity::entryIndex(i, j)];
      } else {
        L[entry_index_ij] = 0.0;
      }
      initImplL<entry_index_ij + 1>(matrix);
    }
  }

  template <int entry_index_ij = 0, class DiagonalBlock>
  [[gnu::always_inline]] void factorImpl(const DiagonalBlock& above) {
    using DiagonalBlockSparsity = typename DiagonalBlock::Sparsity;

    if constexpr (entry_index_ij < nnz) {
      constexpr auto i = Sparsity::entries[entry_index_ij].row_index;
      constexpr auto j = Sparsity::entries[entry_index_ij].col_index;

      static constexpr auto contributions =
          getContributionsMixed<DiagonalBlockSparsity, Sparsity, i, j>();
      auto temp = L[entry_index_ij];
      for (const auto c : contributions) {
        temp -= L[c.entry_index_ik] * above.L[c.entry_index_jk] * above.D[c.k];
      }
      L[entry_index_ij] = temp / above.D[j];

      factorImpl<entry_index_ij + 1>(above);
    }
  }

  template <int entry_index_ij_diag = 0, class DiagonalBlock>
  [[gnu::always_inline]] void contributeImpl(const DiagonalBlock& above,
                                             DiagonalBlock& right) const {
    using DiagonalBlockSparsity = typename DiagonalBlock::Sparsity;
    // TODO: perhaps this should be done in the right DiagonalBlock instead?
    // that way, the whole entry at entry_index_ij_diag would be completed in
    // one go, which would also allow initializing it (with 0 or the original
    // matrix) on the fly

    if constexpr (entry_index_ij_diag < DiagonalBlockSparsity::nnz) {
      constexpr auto i = DiagonalBlockSparsity::entries[entry_index_ij_diag].row_index;
      constexpr auto j = DiagonalBlockSparsity::entries[entry_index_ij_diag].col_index;

      static constexpr auto contributions =
          getContributionsRectangular<Sparsity, i, j>();

      auto temp = right.L[entry_index_ij_diag];
      for (const auto c : contributions) {
        temp -= L[c.entry_index_ik] * L[c.entry_index_jk] * above.D[c.k];
      }
      right.L[entry_index_ij_diag] = temp;

      contributeImpl<entry_index_ij_diag + 1>(above, right);
    }
  }

  template <int entry_index_ij = 0, class DiagonalBlock>
  [[gnu::always_inline]] void contributeImplDiag(const DiagonalBlock& above,
                                                 DiagonalBlock& right) const {
    using DiagonalBlockSparsity = typename DiagonalBlock::Sparsity;
    // TODO: perhaps this should be done in the right DiagonalBlock instead?
    // that way, the whole entry at entry_index_ij_diag would be completed in
    // one go, which would also allow initializing it (with 0 or the original
    // matrix) on the fly

    if constexpr (entry_index_ij < nnz) {
      constexpr auto i = Sparsity::entries[entry_index_ij].row_index;
      constexpr auto j = Sparsity::entries[entry_index_ij].col_index;
      right.D[i] -= square(L[entry_index_ij]) * above.D[j];
      contributeImplDiag<entry_index_ij + 1>(above, right);
    }
  }

  template <int i = 0>
  [[gnu::always_inline]] void forwardSolveImpl(
      const Value* __restrict rhs, Value* __restrict solution) const {
    if constexpr (i < num_rows) {
      auto temp = solution[i];
      constexpr auto row_begin = Sparsity::row_begin_indices[i];
      constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
      for (int entry_index_ij = row_begin; entry_index_ij != row_end;
           ++entry_index_ij) {
        const int j = Sparsity::entries[entry_index_ij].col_index;
        temp -= L[entry_index_ij] * rhs[j];
      }
      solution[i] = temp;
      forwardSolveImpl<i + 1>(rhs, solution);
    }
  }

  template <int i = num_rows - 1>
  [[gnu::always_inline]] void backwardSolveImpl(
      Value* __restrict solution, const Value* __restrict rhs) const {
    if constexpr (i >= 0) {
      const auto temp = rhs[i];
      constexpr auto row_begin = Sparsity::row_begin_indices[i];
      constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
      for (int entry_index_ij = row_begin; entry_index_ij != row_end;
           ++entry_index_ij) {
        const int j = Sparsity::entries[entry_index_ij].col_index;
        solution[j] -= L[entry_index_ij] * temp;
      }
      backwardSolveImpl<i - 1>(solution, rhs);
    }
  }
};

}  // namespace ctldl
