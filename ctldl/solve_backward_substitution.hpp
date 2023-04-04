#pragma once

#include <ctldl/empty_factor_data.hpp>

#include <cstddef>

namespace ctldl {

template <std::size_t i, class FactorData, class Vector>
[[gnu::always_inline]] inline void solveBackwardSubstitutionRow(
    const FactorData& fact, const typename FactorData::Value solution_i,
    Vector& rhs_in_solution_out) {
  using Sparsity = typename FactorData::Sparsity;

  constexpr auto row_begin = std::size_t{Sparsity::row_begin_indices[i]};
  constexpr auto row_end = std::size_t{Sparsity::row_begin_indices[i + 1]};
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = Sparsity::entries[entry_index_ij].col_index;
    rhs_in_solution_out[j] -= fact.L[entry_index_ij] * solution_i;
  }
}

template <std::size_t i, class FactorData, class Vector, class FactorDataLeft,
          class VectorLeft>
[[gnu::always_inline]] inline void solveBackwardSubstitutionImpl(
    const FactorData& diag, Vector& rhs_in_solution_out,
    const FactorDataLeft& left, VectorLeft& rhs_in_solution_out_left) {
  using Sparsity = typename FactorData::Sparsity;
  using SparsityLeft = typename FactorDataLeft::Sparsity;
  static_assert(Sparsity::num_rows == SparsityLeft::num_rows);

  const auto solution_i = rhs_in_solution_out[i];
  solveBackwardSubstitutionRow<i>(diag, solution_i, rhs_in_solution_out);
  solveBackwardSubstitutionRow<i>(left, solution_i, rhs_in_solution_out_left);
  if constexpr (i > 0) {
    solveBackwardSubstitutionImpl<i - 1>(diag, rhs_in_solution_out, left,
                                         rhs_in_solution_out_left);
  }
}

template <class FactorData, class Vector, class FactorDataLeft,
          class VectorLeft>
void solveBackwardSubstitution(const FactorData& diag,
                               Vector& rhs_in_solution_out,
                               const FactorDataLeft& left,
                               VectorLeft& rhs_in_solution_out_left) {
  constexpr auto num_rows = std::size_t{FactorData::Sparsity::num_rows};
  if constexpr (num_rows > 0) {
    solveBackwardSubstitutionImpl<num_rows - 1>(diag, rhs_in_solution_out, left,
                                                rhs_in_solution_out_left);
  }
}

template <class FactorData, class Vector>
void solveBackwardSubstitution(const FactorData& diag,
                               Vector& rhs_in_solution_out) {
  using Value = typename FactorData::Value;

  const EmptyFactorData<FactorData::Sparsity::num_rows, Value> empty_left;
  std::array<Value, 0> empty_rhs_in_solution_out_left;
  solveBackwardSubstitution(diag, rhs_in_solution_out, empty_left,
                            empty_rhs_in_solution_out_left);
}

}  // namespace ctldl
