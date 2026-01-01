#pragma once

#include <ctldl/factor_data/empty_factor_data_left.hpp>
#include <ctldl/utility/unroll.hpp>

#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t i, class FactorData, class Vector, class ValueOut>
[[gnu::always_inline]] inline auto solveForwardSubstitutionRow(
    const FactorData& fact, const Vector& solution, const ValueOut init) {
  constexpr auto& sparsity = FactorData::sparsity;

  constexpr auto row_begin = std::size_t{sparsity.rowBeginIndices()[i]};
  constexpr auto row_end = std::size_t{sparsity.rowBeginIndices()[i + 1]};
  auto value = init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = sparsity.entries()[entry_index_ij].col_index;
    const auto j_orig = FactorData::origColIndex(j);
    value -= static_cast<ValueOut>(fact.L[entry_index_ij]) *
             static_cast<ValueOut>(solution[j_orig]);
  }
  return value;
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
  constexpr auto num_rows = std::size_t{FactorData::sparsity.numRows()};
  unroll<0, num_rows>([&](const auto i) {
    constexpr auto i_orig = FactorData::origRowIndex(i);
    partial_solution[i_orig] = solveForwardSubstitutionRow<i>(
        factor_block, solution, partial_solution[i_orig]);
  });
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
  static_assert(FactorData::sparsity.numRows() ==
                FactorDataLeft::sparsity.numRows());
  static_assert(FactorData::permutation == FactorDataLeft::permutation_row);

  constexpr auto num_rows = std::size_t{FactorData::sparsity.numRows()};
  unroll<0, num_rows>([&](const auto i) {
    constexpr auto i_orig = FactorData::origRowIndex(i);
    auto solution_i = rhs_in_solution_out[i_orig];
    solution_i = solveForwardSubstitutionRow<i>(left, solution_left, solution_i);
    solution_i = solveForwardSubstitutionRow<i>(diag, rhs_in_solution_out, solution_i);
    rhs_in_solution_out[i_orig] = solution_i;
  });
}

template <class FactorData, class Vector>
void solveForwardSubstitution(const FactorData& diag,
                              Vector&& rhs_in_solution_out) {
  solveForwardSubstitution(diag, rhs_in_solution_out, rhs_in_solution_out);
}

}  // namespace ctldl
