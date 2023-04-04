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

  template <class Matrix, class DiagonalBlock>
  void factor(const Matrix& matrix, const DiagonalBlock& above) {
    factorImpl(matrix, above);
  }

  template <class Rhs>
  void forwardSolve(const Rhs& rhs, Rhs& solution) const {
    forwardSolveImpl(rhs.data(), solution.data());
  }

  template <class Rhs>
  void backwardSolve(Rhs& solution, const Rhs& rhs) const {
    backwardSolveImpl(solution.data(), rhs.data());
  }

  std::array<Value, Sparsity::nnz> L;

 private:
  template <int entry_index_ij = 0, class Matrix, class DiagonalBlock>
  [[gnu::always_inline]] void factorImpl(const Matrix& matrix,
                                         const DiagonalBlock& above) {
    using DiagonalBlockSparsity = typename DiagonalBlock::Sparsity;

    if constexpr (entry_index_ij < nnz) {
      constexpr auto i = Sparsity::entries[entry_index_ij].row_index;
      constexpr auto j = Sparsity::entries[entry_index_ij].col_index;

      static constexpr auto contributions =
          getContributionsMixed<DiagonalBlockSparsity, Sparsity, i, j>();
      auto Lij = matrix.get(i, j);
      for (const auto c : contributions) {
        Lij -= L[c.entry_index_ik] * above.L[c.entry_index_jk] * above.D[c.k];
      }
      L[entry_index_ij] = Lij / above.D[j];

      factorImpl<entry_index_ij + 1>(matrix, above);
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
