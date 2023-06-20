#pragma once

#include <ctldl/empty_factor_data_left.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

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
    const auto j_orig = FactorData::permutation_col[j];
    rhs_in_solution_out[j_orig] -= fact.L[entry_index_ij] * solution_i;
  }
}

template <std::size_t i, class Vector, class FactorDataLeft, class VectorLeft>
[[gnu::always_inline]] inline void solveBackwardSubstitutionImpl(
    Vector& rhs_in_solution_out, const FactorDataLeft& left,
    const VectorLeft& rhs_in_solution_out_left) {
  using Value = typename FactorDataLeft::Value;

  constexpr auto i_orig_left = FactorDataLeft::permutation_row[i];
  const auto solution_i_left =
      static_cast<Value>(rhs_in_solution_out_left[i_orig_left]);
  solveBackwardSubstitutionRow<i>(left, solution_i_left, rhs_in_solution_out);
}

template <std::size_t i, class FactorData, class Vector>
[[gnu::always_inline]] inline void solveBackwardSubstitutionImpl(
    const FactorData& diag, Vector& rhs_in_solution_out) {
  using Value = typename FactorData::Value;

  constexpr auto i_orig = FactorData::permutation[i];
  const auto solution_i = static_cast<Value>(rhs_in_solution_out[i_orig]);
  solveBackwardSubstitutionRow<i>(diag, solution_i, rhs_in_solution_out);
}

template <std::size_t... RowIndices, class Vector, class FactorDataLeft,
          class VectorLeft>
void solveBackwardSubstitutionImpl(Vector& rhs_in_solution_out,
                                   const FactorDataLeft& left,
                                   const VectorLeft& rhs_in_solution_out_left,
                                   std::index_sequence<RowIndices...>) {
  (solveBackwardSubstitutionImpl<RowIndices>(rhs_in_solution_out, left,
                                             rhs_in_solution_out_left),
   ...);
}

template <std::size_t... RowIndices, class FactorData, class Vector>
void solveBackwardSubstitutionImpl(const FactorData& diag,
                                   Vector& rhs_in_solution_out,
                                   std::index_sequence<RowIndices...>) {
  (solveBackwardSubstitutionImpl<RowIndices>(diag, rhs_in_solution_out), ...);
}

template <class FactorData, class Vector, class FactorDataLeft,
          class VectorLeft>
void solveBackwardSubstitution(const FactorData& diag,
                               Vector& rhs_in_solution_out,
                               const FactorDataLeft& left,
                               const VectorLeft& rhs_in_solution_out_left) {
  constexpr auto num_rows = std::size_t{FactorData::Sparsity::num_rows};
  solveBackwardSubstitutionImpl(rhs_in_solution_out, left,
                                rhs_in_solution_out_left,
                                makeIndexSequenceReversed<0, num_rows>());
  solveBackwardSubstitutionImpl(diag, rhs_in_solution_out,
                                makeIndexSequenceReversed<0, num_rows>());
}

template <class FactorData, class Vector>
void solveBackwardSubstitution(const FactorData& diag,
                               Vector& rhs_in_solution_out) {
  using Value = typename FactorData::Value;

  constexpr EmptyFactorDataLeft<FactorData> empty_left;
  std::array<Value, 0> empty_rhs_in_solution_out_left;
  solveBackwardSubstitution(diag, rhs_in_solution_out, empty_left,
                            empty_rhs_in_solution_out_left);
}

}  // namespace ctldl
