#pragma once

#include <ctldl/empty_factor_data_left.hpp>

#include <cstddef>

namespace ctldl {

template <std::size_t i, class FactorData, class Vector>
[[gnu::always_inline]] inline auto solveForwardSubstitutionRow(
    const FactorData& fact, const Vector& solution,
    const typename FactorData::Value init) {
  using Sparsity = typename FactorData::Sparsity;

  constexpr auto row_begin = std::size_t{Sparsity::row_begin_indices[i]};
  constexpr auto row_end = std::size_t{Sparsity::row_begin_indices[i + 1]};
  auto value = init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = Sparsity::entries[entry_index_ij].col_index;
    value -= fact.L[entry_index_ij] * solution[j];
  }
  return value;
}

template <std::size_t i = 0, class FactorData, class Vector,
          class FactorDataLeft, class VectorLeft>
[[gnu::always_inline]] inline void solveForwardSubstitutionImpl(
    const FactorData& diag, Vector& rhs_in_solution_out,
    const FactorDataLeft& left, const VectorLeft& solution_left) {
  using Sparsity = typename FactorData::Sparsity;
  using SparsityLeft = typename FactorDataLeft::Sparsity;
  static_assert(Sparsity::num_rows == SparsityLeft::num_rows);

  if constexpr (i < Sparsity::num_rows) {
    auto temp = rhs_in_solution_out[i];
    temp = solveForwardSubstitutionRow<i>(left, solution_left, temp);
    temp = solveForwardSubstitutionRow<i>(diag, rhs_in_solution_out, temp);
    rhs_in_solution_out[i] = temp;
    solveForwardSubstitutionImpl<i + 1>(diag, rhs_in_solution_out, left,
                                        solution_left);
  }
}

template <class FactorData, class Vector, class FactorDataLeft,
          class VectorLeft>
void solveForwardSubstitution(const FactorData& diag,
                              Vector& rhs_in_solution_out,
                              const FactorDataLeft& left,
                              const VectorLeft& solution_left) {
  solveForwardSubstitutionImpl(diag, rhs_in_solution_out, left, solution_left);
}

template <class FactorData, class Vector>
void solveForwardSubstitution(const FactorData& diag,
                              Vector& rhs_in_solution_out) {
  using Value = typename FactorData::Value;

  constexpr EmptyFactorDataLeft<FactorData> empty_left;
  constexpr std::array<Value, 0> empty_solution_left;
  solveForwardSubstitution(diag, rhs_in_solution_out, empty_left,
                           empty_solution_left);
}

}  // namespace ctldl
