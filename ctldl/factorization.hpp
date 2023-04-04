#pragma once

#include <ctldl/factorize.hpp>
#include <ctldl/sparsity/filled_in_sparsity.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>
#include <ctldl/utility/square.hpp>

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

template <class OriginalSparsity, class Value_>
class Factorization {
 public:
  static_assert(OriginalSparsity::num_rows == OriginalSparsity::num_cols);
  static constexpr auto dim = std::size_t{OriginalSparsity::num_rows};
  using Sparsity = SparsityCSR<FilledInSparsity<OriginalSparsity>>;
  static constexpr auto nnz = std::size_t{Sparsity::nnz};
  using Value = Value_;

  template <class Matrix>
  void factor(const Matrix& matrix) {
    factorize(*this, matrix);
  }

  template <class Matrix, class FactorDataLeft>
  void factor(const Matrix& matrix, const FactorDataLeft& left) {
    factorize(*this, matrix, left);
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
  template <std::size_t i = 0>
  [[gnu::always_inline]] void forwardSolveImpl(
      Value* __restrict rhs_in_solution_out) const {
    if constexpr (i < dim) {
      auto temp = rhs_in_solution_out[i];
      constexpr auto row_begin = std::size_t{Sparsity::row_begin_indices[i]};
      constexpr auto row_end = std::size_t{Sparsity::row_begin_indices[i + 1]};
      for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
           ++entry_index_ij) {
        const std::size_t j = Sparsity::entries[entry_index_ij].col_index;
        temp -= L[entry_index_ij] * rhs_in_solution_out[j];
      }
      rhs_in_solution_out[i] = temp;
      forwardSolveImpl<i + 1>(rhs_in_solution_out);
    }
  }

  void diagonalSolveImpl(Value* __restrict rhs_in_solution_out) const {
    for (std::size_t i = 0; i < dim; ++i) {
      rhs_in_solution_out[i] /= D[i];
    }
  }

  template <std::size_t i = dim>
  [[gnu::always_inline]] void backwardSolveImpl(
      Value* __restrict rhs_in_solution_out) const {
    if constexpr (i > 0) {
      const auto temp = rhs_in_solution_out[i - 1];
      constexpr auto row_begin = std::size_t{Sparsity::row_begin_indices[i - 1]};
      constexpr auto row_end = std::size_t{Sparsity::row_begin_indices[i]};
      for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
           ++entry_index_ij) {
        const std::size_t j = Sparsity::entries[entry_index_ij].col_index;
        rhs_in_solution_out[j] -= L[entry_index_ij] * temp;
      }
      backwardSolveImpl<i - 1>(rhs_in_solution_out);
    }
  }
};

}  // namespace ctldl
