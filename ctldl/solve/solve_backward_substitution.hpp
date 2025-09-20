#pragma once

#include <ctldl/factor_data/empty_factor_data_left.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

#include <cstddef>

namespace ctldl {

template <std::size_t i, class FactorData, class ValueSolution, class Vector>
[[gnu::always_inline]] inline void solveBackwardSubstitutionRow(
    const FactorData& fact, const ValueSolution solution_i,
    Vector& partial_solution) {
  using ValueOut = std::remove_reference_t<decltype(partial_solution[0])>;
  constexpr auto& sparsity = FactorData::sparsity;

  constexpr auto row_begin = std::size_t{sparsity.row_begin_indices[i]};
  constexpr auto row_end = std::size_t{sparsity.row_begin_indices[i + 1]};
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = sparsity.entries[entry_index_ij].col_index;
    const auto j_orig = FactorData::origColIndex(j);
    partial_solution[j_orig] -= static_cast<ValueOut>(fact.L[entry_index_ij]) *
                                static_cast<ValueOut>(solution_i);
  }
}

template <std::size_t i, class FactorData, class VectorSolution,
          class VectorPartialSolution>
[[gnu::always_inline]] inline void solveBackwardSubstitutionImpl(
    const FactorData& factor_block, const VectorSolution& solution,
    VectorPartialSolution&& partial_solution) {
  constexpr auto i_orig = FactorData::origRowIndex(i);
  solveBackwardSubstitutionRow<i>(factor_block, solution[i_orig],
                                  partial_solution);
}

template <std::size_t... RowIndices, class FactorData, class VectorSolution,
          class VectorPartialSolution>
void solveBackwardSubstitutionImpl(const FactorData& factor_block,
                                   const VectorSolution& solution,
                                   VectorPartialSolution&& partial_solution,
                                   std::index_sequence<RowIndices...>) {
  (solveBackwardSubstitutionImpl<RowIndices>(factor_block, solution,
                                             partial_solution),
   ...);
}

/**
 * Performs a single block-operation of the backward substitution with the
 * factor of a block-matrix.
 *
 * For the example
 * \code
 * [*  *  *  *  *]  [*]   [*]
 * [   *  * F^T *]  [*]   [r]
 * [      *  *  *]  [*] = [*]
 * [         *  *]  [s]   [*]
 * [            *]  [*]   [*]
 * \endcode
 *
 * one of these block operations would be subtracting the product of F^T
 * and s from r. That is exactly what this function does with
 * F: \p factor_block, s: \p solution and r: \p partial_solution. So it
 * basically performs the matrix-vector operation
 * partial_solution -= transposed(factor_block) * solution.
 *
 * The underlying permutation of the original matrix that was factorized is
 * applied internally, so \p partial_solution and \p solution should be given
 * unpermuted, corresponding to the original matrix order.
 */
template <class FactorData, class VectorSolution, class VectorPartialSolution>
void solveBackwardSubstitution(const FactorData& factor_block,
                               const VectorSolution& solution,
                               VectorPartialSolution&& partial_solution) {
  constexpr auto num_rows = std::size_t{FactorData::sparsity.num_rows};
  solveBackwardSubstitutionImpl(factor_block, solution, partial_solution,
                                makeIndexSequenceReversed<0, num_rows>());
}

/**
 * Performs backward substitution for a single diagonal-block of the factor of
 * a block-matrix.
 *
 * This is just a special overload for the block being on the diagonal of the
 * factor. For this case, the values that are computed in \p partial_solution
 * in one entry become the final solution and are then needed soon after for
 * computing an earlier entry. Therefore it makes sense to do the operation
 * in-place. So \p partial_solution should contain the partial solution of the
 * backward substitution done with all non-diagonal blocks, such that it returns
 * the complete solution on return.
 */
template <class FactorData, class Vector>
void solveBackwardSubstitution(const FactorData& factor_block_diagonal,
                               Vector&& partial_solution) {
  solveBackwardSubstitution(factor_block_diagonal, partial_solution,
                            partial_solution);
}

}  // namespace ctldl
