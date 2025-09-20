#pragma once

#include <ctldl/factor_data/empty_factor_data_left.hpp>

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t i, class FactorData, class Vector, class ValueOut>
[[gnu::always_inline]] inline auto solveForwardSubstitutionRow(
    const FactorData& fact, const Vector& solution, const ValueOut init) {
  constexpr auto& sparsity = FactorData::sparsity;

  constexpr auto row_begin = std::size_t{sparsity.row_begin_indices[i]};
  constexpr auto row_end = std::size_t{sparsity.row_begin_indices[i + 1]};
  auto value = init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = sparsity.entries[entry_index_ij].col_index;
    const auto j_orig = FactorData::origColIndex(j);
    value -= static_cast<ValueOut>(fact.L[entry_index_ij]) *
             static_cast<ValueOut>(solution[j_orig]);
  }
  return value;
}

template <std::size_t i, class FactorData, class VectorSolution,
          class VectorPartialSolution>
[[gnu::always_inline]] inline auto solveForwardSubstitutionImpl(
    const FactorData& factor_block, const VectorSolution& solution,
    const VectorPartialSolution& partial_solution) {
  constexpr auto i_orig = FactorData::origRowIndex(i);
  return solveForwardSubstitutionRow<i>(factor_block, solution,
                                        partial_solution[i_orig]);
}

template <std::size_t... RowIndices, class FactorData, class VectorSolution,
          class VectorPartialSolution>
void solveForwardSubstitutionImpl(const FactorData& factor_block,
                                  const VectorSolution& solution,
                                  VectorPartialSolution&& partial_solution,
                                  std::index_sequence<RowIndices...>) {
  ((partial_solution[FactorData::origRowIndex(RowIndices)] =
        solveForwardSubstitutionImpl<RowIndices>(factor_block, solution,
                                                 partial_solution)),
   ...);
}

/**
 * Performs a single block-operation of the forward substitution with the
 * factor of a block-matrix.
 *
 * For the example
 * \code
 * [*            ]  [*]   [*]
 * [*  *         ]  [s]   [*]
 * [*  *  *      ]  [*] = [*]
 * [*  F  *  *   ]  [*]   [r]
 * [*  *  *  *  *]  [*]   [*]
 * \endcode
 *
 * one of these block operations would be subtracting the product of F and s
 * from r. That is exactly what this function does with F: \p factor_block,
 * s: \p solution and r: \p partial_solution. So itbasically performs the
 * matrix-vector operation partial_solution -= factor_block * solution.
 *
 * The underlying permutation of the original matrix that was factorized is
 * applied internally, so \p partial_solution and \p solution should be given
 * unpermuted, corresponding to the original matrix order.
 */
template <class FactorData, class VectorSolution, class VectorPartialSolution>
void solveForwardSubstitution(const FactorData& factor_block,
                              const VectorSolution& solution,
                              VectorPartialSolution&& partial_solution) {
  constexpr auto num_rows = std::size_t{FactorData::sparsity.num_rows};
  solveForwardSubstitutionImpl(factor_block, solution, partial_solution,
                               std::make_index_sequence<num_rows>());
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

  constexpr auto i_orig = FactorData::origRowIndex(i);
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
  ((rhs_in_solution_out[FactorData::origRowIndex(RowIndices)] =
        solveForwardSubstitutionImpl<RowIndices>(diag, rhs_in_solution_out,
                                                 left, solution_left)),
   ...);
}

/**
 * Optimized overload for forward substitution of a block-row with exactly one
 * off-diagonal and one diagonal block in one go.
 *
 * For the example
 * \code
 * [*            ]  [*]   [*]
 * [*  *         ]  [s]   [*]
 * [*  *  *      ]  [*] = [*]
 * [0  F  0  D   ]  [*]   [r]
 * [*  *  *  *  *]  [*]   [*]
 * \endcode
 *
 * this combined block operation would be subtracting the product of F and s and
 * the product of D and s from r. That is exactly what this function does with
 * F: \p left, D: \p diag, s: \p solution_left and r: \p rhs_in_solution_out. So
 * it basically performs the matrix-vector operation
 * rhs_in_solution_out -= left * solution_left + diag * rhs_in_solution_out.
 *
 * The underlying permutation of the original matrix that was factorized is
 * applied internally, so \p solution_left and \p rhs_in_solution_out should be
 * given unpermuted, corresponding to the original matrix order.
 */
template <class FactorData, class Vector, class FactorDataLeft,
          class VectorLeft>
void solveForwardSubstitution(const FactorData& diag,
                              Vector&& rhs_in_solution_out,
                              const FactorDataLeft& left,
                              const VectorLeft& solution_left) {
  constexpr auto num_rows = std::size_t{FactorData::sparsity.num_rows};
  solveForwardSubstitutionImpl(diag, rhs_in_solution_out, left, solution_left,
                               std::make_index_sequence<num_rows>());
}

template <class FactorData, class Vector>
void solveForwardSubstitution(const FactorData& diag,
                              Vector&& rhs_in_solution_out) {
  solveForwardSubstitution(diag, rhs_in_solution_out, rhs_in_solution_out);
}

}  // namespace ctldl
