#pragma once

#include <ctldl/sparsity/filled_in_sparsity.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/utility/square.hpp>

#include <algorithm>
#include <array>

namespace ctldl {

template <class OriginalSparsity, class Value_>
class Factorization {
 public:
  static_assert(OriginalSparsity::num_rows == OriginalSparsity::num_cols);
  static constexpr int dim = OriginalSparsity::num_rows;
  using Sparsity = SparsityCSR<FilledInSparsity<OriginalSparsity>>;
  static constexpr int nnz = Sparsity::nnz;
  using Value = Value_;

  template <class Matrix>
  void init(const Matrix& matrix) {
    initImplD(matrix);
    initImplL(matrix);
  }

  void factor() {
    factorImpl();
  }

  template <class Rhs>
  void forwardSolve(Rhs& rhs_in_solution_out) const {
    forwardSolveImpl(rhs_in_solution_out.data());
  }

  template <class Rhs>
  void diagonalSolve(Rhs& rhs_in_solution_out) const {
    diagonalSolveImpl(rhs_in_solution_out.data());
  }

  template <class Rhs>
  void backwardSolve(Rhs& rhs_in_solution_out) const {
    backwardSolveImpl(rhs_in_solution_out.data());
  }

  std::array<Value, nnz> L;
  std::array<Value, dim> D;

 private:
  template <int i = 0, class Matrix>
  [[gnu::always_inline]] void initImplD(const Matrix& matrix) {
    if constexpr (i < dim) {
      if constexpr (Matrix::Sparsity::is_nonzero[i][i]) {
        D[i] = matrix.m_values[Matrix::Sparsity::entryIndex(i, i)];
      } else {
        D[i] = 0.0;
      }
      initImplD<i + 1>(matrix);
    }
  }

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

  template <int i = 0>
  [[gnu::always_inline]] void factorImpl() {
    if constexpr (i < dim) {
      D[i] = factorImplRow<i>(D[i]);
      factorImpl<i + 1>();
    }
  }

  template <int i, int entry_index_ij = Sparsity::row_begin_indices[i]>
  [[gnu::always_inline]] Value factorImplRow(const Value Di_prev) {
    constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
    if constexpr (entry_index_ij < row_end) {
      constexpr auto j = Sparsity::entries[entry_index_ij].col_index;

      static constexpr auto contributions =
          getContributionsLowerTriangle<Sparsity, i, j>();
      auto temp = L[entry_index_ij];
      for (const auto c : contributions) {
        temp -= L[c.entry_index_ik] * L[c.entry_index_jk] * D[c.k];
      }
      const auto Lij = temp / D[j];
      L[entry_index_ij] = Lij;

      const auto Di = factorImplRow<i, entry_index_ij + 1>(Di_prev);
      return Di - square(Lij) * D[j];
    }
    return Di_prev;
  }

  template <int i = 0>
  [[gnu::always_inline]] void forwardSolveImpl(
      Value* __restrict rhs_in_solution_out) const {
    if constexpr (i < dim) {
      auto temp = rhs_in_solution_out[i];
      constexpr auto row_begin = Sparsity::row_begin_indices[i];
      constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
      for (int entry_index_ij = row_begin; entry_index_ij != row_end;
           ++entry_index_ij) {
        const int j = Sparsity::entries[entry_index_ij].col_index;
        temp -= L[entry_index_ij] * rhs_in_solution_out[j];
      }
      rhs_in_solution_out[i] = temp;
      forwardSolveImpl<i + 1>(rhs_in_solution_out);
    }
  }

  void diagonalSolveImpl(Value* __restrict rhs_in_solution_out) const {
    for (int i = 0; i < dim; ++i) {
      rhs_in_solution_out[i] /= D[i];
    }
  }

  template <int i = dim - 1>
  [[gnu::always_inline]] void backwardSolveImpl(
      Value* __restrict rhs_in_solution_out) const {
    if constexpr (i >= 0) {
      const auto temp = rhs_in_solution_out[i];
      constexpr auto row_begin = Sparsity::row_begin_indices[i];
      constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
      for (int entry_index_ij = row_begin; entry_index_ij != row_end;
           ++entry_index_ij) {
        const int j = Sparsity::entries[entry_index_ij].col_index;
        rhs_in_solution_out[j] -= L[entry_index_ij] * temp;
      }
      backwardSolveImpl<i - 1>(rhs_in_solution_out);
    }
  }
};

}  // namespace ctldl
