#pragma once

#include <ctldl/factor_data/empty_factor_data_left.hpp>

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t i, class FactorData, class Vector>
[[gnu::always_inline]] inline auto solveForwardSubstitutionRow(
    const FactorData& fact, const Vector& solution,
    const typename FactorData::Value init) {
  constexpr auto& sparsity = FactorData::sparsity;
  using Value = typename FactorData::Value;

  constexpr auto row_begin = std::size_t{sparsity.row_begin_indices[i]};
  constexpr auto row_end = std::size_t{sparsity.row_begin_indices[i + 1]};
  auto value = init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = sparsity.entries[entry_index_ij].col_index;
    const auto j_orig = FactorData::permutation_col[j];
    value -= fact.L[entry_index_ij] * static_cast<Value>(solution[j_orig]);
  }
  return value;
}

template <std::size_t i, class FactorData, class Vector, class FactorDataLeft,
          class VectorLeft>
[[gnu::always_inline]] inline auto solveForwardSubstitutionImpl(
    const FactorData& diag, const Vector& rhs_in_solution_out,
    const FactorDataLeft& left, const VectorLeft& solution_left) {
  static_assert(FactorData::sparsity.num_rows ==
                FactorDataLeft::sparsity.num_rows);
  static_assert(FactorData::permutation == FactorDataLeft::permutation_row);
  using Value = typename FactorData::Value;

  constexpr auto i_orig = FactorData::permutation[i];
  auto solution_i = static_cast<Value>(rhs_in_solution_out[i_orig]);
  solution_i = solveForwardSubstitutionRow<i>(left, solution_left, solution_i);
  solution_i = solveForwardSubstitutionRow<i>(diag, rhs_in_solution_out, solution_i);
  return solution_i;
}

template <std::size_t... RowIndices, class FactorData, class Vector,
          class FactorDataLeft, class VectorLeft>
void solveForwardSubstitutionImpl(const FactorData& diag,
                                  Vector& rhs_in_solution_out,
                                  const FactorDataLeft& left,
                                  const VectorLeft& solution_left,
                                  std::index_sequence<RowIndices...>) {
  ((rhs_in_solution_out[FactorData::permutation[RowIndices]] =
        solveForwardSubstitutionImpl<RowIndices>(diag, rhs_in_solution_out,
                                                 left, solution_left)),
   ...);
}

template <class FactorData, class Vector, class FactorDataLeft,
          class VectorLeft>
void solveForwardSubstitution(const FactorData& diag,
                              Vector& rhs_in_solution_out,
                              const FactorDataLeft& left,
                              const VectorLeft& solution_left) {
  constexpr auto num_rows = std::size_t{FactorData::sparsity.num_rows};
  solveForwardSubstitutionImpl(diag, rhs_in_solution_out, left, solution_left,
                               std::make_index_sequence<num_rows>());
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
