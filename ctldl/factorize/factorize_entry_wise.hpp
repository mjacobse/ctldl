#pragma once

#include <ctldl/factor_data/empty_factor_data_diagonal.hpp>
#include <ctldl/factor_data/empty_factor_data_left.hpp>
#include <ctldl/factorize/factor_initialization.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/factorize/regularization.hpp>
#include <ctldl/matrix/empty_matrix_input.hpp>
#include <ctldl/permutation/permuted_entry_lower_triangle.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/is_sparsity_subset.hpp>
#include <ctldl/symbolic/is_chordal_blocked.hpp>
#include <ctldl/utility/make_index_sequence.hpp>
#include <ctldl/utility/square.hpp>

#include <array>
#include <cstddef>
#include <utility>

namespace ctldl {

/**
 * Factorize the entry given by \p entry_index in block \p factor21 (L21) using
 * the initial values for that block given by \p init21 and the already
 * finished block \p factor11 (L11).
 *
 * [*            ]
 * [* L11        ]
 * [*  *  *      ]
 * [* L21 *  *   ]
 * [*  *  *  *  *]
 */
template <std::size_t entry_index_ij, class FactorData21, class Init21,
          class FactorData11>
[[gnu::always_inline]] inline void factorEntryWiseSubdiagonalImplRow(
    FactorData21& factor21, const Init21& init21, const FactorData11& factor11) {
  constexpr auto& sparsity11 = FactorData11::sparsity;
  constexpr auto& sparsity21 = FactorData21::sparsity;

  constexpr auto i = std::size_t{sparsity21.entries()[entry_index_ij].row_index};
  constexpr auto j = std::size_t{sparsity21.entries()[entry_index_ij].col_index};

  auto Lij = getInitialFactorValueL<entry_index_ij>(init21, factor21);

  static constexpr auto contributions =
      getContributionsMixed<sparsity11, sparsity21, i, j>();
  for (const auto c : contributions) {
    Lij -= factor21.L[c.entry_index_ik] * factor11.L[c.entry_index_jk] *
           factor11.D[c.k];
  }
  factor21.L[entry_index_ij] = Lij / factor11.D[j];
}

template <std::size_t... EntryIndices, class FactorData21, class Init21,
          class FactorData11>
[[gnu::always_inline]] inline void factorEntryWiseSubdiagonalImplRow(
    FactorData21& factor21, const Init21& init21, const FactorData11& factor11,
    std::index_sequence<EntryIndices...>) {
  (factorEntryWiseSubdiagonalImplRow<EntryIndices>(factor21, init21, factor11),
   ...);
}

template <std::size_t i, class FactorData21, class Init21, class FactorData11>
[[gnu::always_inline]] inline void factorizeEntryWiseSubdiagonalImpl(
    FactorData21& factor21, const Init21& init21, const FactorData11& factor11) {
  constexpr auto& sparsity21 = FactorData21::sparsity;
  constexpr auto row_begin = std::size_t{sparsity21.rowBeginIndices()[i]};
  constexpr auto row_end = std::size_t{sparsity21.rowBeginIndices()[i + 1]};
  factorEntryWiseSubdiagonalImplRow(factor21, init21, factor11,
                                    makeIndexSequence<row_begin, row_end>());
}

template <std::size_t... RowIndices, class FactorData21, class Init21,
          class FactorData11>
void factorizeEntryWiseSubdiagonalImpl(FactorData21& factor21,
                                       const Init21& init21,
                                       const FactorData11& factor11,
                                       std::index_sequence<RowIndices...>) {
  (factorizeEntryWiseSubdiagonalImpl<RowIndices>(factor21, init21, factor11),
   ...);
}

template <std::size_t i, class FactorData, class FactorDataDiag>
[[gnu::always_inline]] inline auto applyContributionsRowDiagonal(
    const FactorData& fact, const FactorDataDiag& diag,
    const typename FactorData::Value value_init) {
  constexpr auto& sparsity = FactorData::sparsity;

  constexpr auto row_begin = std::size_t{sparsity.rowBeginIndices()[i]};
  constexpr auto row_end = std::size_t{sparsity.rowBeginIndices()[i + 1]};
  auto value = value_init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = sparsity.entries()[entry_index_ij].col_index;
    value -= square(fact.L[entry_index_ij]) * diag.D[j];
  }
  return value;
}

template <class FactorData, class FactorDataDiag, std::size_t num_contributions>
[[gnu::always_inline]] inline auto applyContributions(
    const FactorData& fact, const FactorDataDiag& diag,
    const std::array<Contribution, num_contributions>& contributions,
    const typename FactorData::Value value_init) {
  auto value = value_init;
  for (const auto c : contributions) {
    value -= fact.L[c.entry_index_ik] * fact.L[c.entry_index_jk] * diag.D[c.k];
  }
  return value;
}

/**
 * Factorize the entry given by \p entry_index in block \p factor22 (L22) using
 * the initial values for that block given by \p init22 and the already finished
 * blocks \p factor11 (L11) and \p factor21 (L21).
 *
 * [*            ]
 * [* L11        ]
 * [*  *  *      ]
 * [* L21 * L22  ]
 * [*  *  *  *  *]
 *
 * \return The contribution of the finished entry to the diagonal entry in its
 * row.
 */
template <std::size_t entry_index_ij, class FactorData22, class Init22,
          class FactorData21, class FactorData11>
[[gnu::always_inline]] inline auto factorizeEntryWiseImplRow(
    FactorData22& factor22, const Init22& init22, const FactorData21& factor21,
    const FactorData11& factor11) {
  constexpr auto& sparsity21 = FactorData21::sparsity;
  constexpr auto& sparsity22 = FactorData22::sparsity;

  constexpr auto i = std::size_t{sparsity22.entries()[entry_index_ij].row_index};
  constexpr auto j = std::size_t{sparsity22.entries()[entry_index_ij].col_index};

  auto Lij = getInitialFactorValueL<entry_index_ij>(init22, factor22);

  static constexpr auto contributions21 =
      getContributionsRectangular<sparsity21, i, j>();
  static constexpr auto contributions22 =
      getContributionsLowerTriangle<sparsity22, i, j>();

  Lij = applyContributions(factor21, factor11, contributions21, Lij);
  Lij = applyContributions(factor22, factor22, contributions22, Lij);
  Lij /= factor22.D[j];
  factor22.L[entry_index_ij] = Lij;

  return square(Lij) * factor22.D[j];
}

template <std::size_t... EntryIndices, class FactorData22, class Init22,
          class FactorData21, class FactorData11>
[[gnu::always_inline]] inline auto factorizeEntryWiseImplRow(
    FactorData22& factor22, const Init22& init22, const FactorData21& factor21,
    const FactorData11& factor11, const typename FactorData22::Value Di_init,
    std::index_sequence<EntryIndices...>) {
  return (Di_init - ... -
          factorizeEntryWiseImplRow<EntryIndices>(factor22, init22, factor21,
                                                  factor11));
}

template <std::size_t i, class FactorData22, class Init22, class FactorData21,
          class FactorData11>
[[gnu::always_inline]] inline void factorizeEntryWiseImpl(
    FactorData22& factor22, const Init22& init22, const FactorData21& factor21,
    const FactorData11& factor11, const Regularization auto& regularization) {
  constexpr auto& sparsity21 = FactorData21::sparsity;
  constexpr auto& sparsity22 = FactorData22::sparsity;
  static_assert(sparsity22.numRows() == sparsity21.numRows());

  constexpr auto row_begin = std::size_t{sparsity22.rowBeginIndices()[i]};
  constexpr auto row_end = std::size_t{sparsity22.rowBeginIndices()[i + 1]};

  auto Di = getInitialFactorValueD<i>(init22, factor22);
  Di = applyContributionsRowDiagonal<i>(factor21, factor11, Di);
  Di = factorizeEntryWiseImplRow(factor22, init22, factor21, factor11, Di,
                                 makeIndexSequence<row_begin, row_end>());
  Di = regularization.regularize(Di, i);
  factor22.D[i] = Di;
}

template <std::size_t... RowIndices, class FactorData22, class Init22,
          class FactorData21, class FactorData11>
void factorizeEntryWiseImpl(FactorData22& factor22, const Init22& init22,
                            const FactorData21& factor21,
                            const FactorData11& factor11,
                            const Regularization auto& regularization,
                            std::index_sequence<RowIndices...>) {
  (factorizeEntryWiseImpl<RowIndices>(factor22, init22, factor21, factor11,
                                      regularization),
   ...);
}

/**
 * Factorize blocks \p factor21 (L21) and \p factor22 (L22) using inital values
 * given by \p init21 and \p init22 for those blocks and the already finished
 * block \p factor11 (L11) in an entry-wise way.
 *
 * [*            ]
 * [* L11        ]
 * [*  *  *      ]
 * [* L21 * L22  ]
 * [*  *  *  *  *]
 *
 * Note that this implicitly assumes that no other blocks have influence on L21
 * and L22, since no partial result can be returned that other contributions
 * could be added to. Instead the final result is computed in one go directly
 * starting with the input values.
 */
template <class FactorData11, class Init21, class Init22, class FactorData21,
          class FactorData22>
void factorizeEntryWise(const FactorData11& factor11, const Init21& init21,
                        const Init22& init22, FactorData21& factor21,
                        FactorData22& factor22,
                        const Regularization auto& regularization) {
  static_assert(isChordalBlocked(FactorData11::sparsity, FactorData21::sparsity,
                                 FactorData22::sparsity));
  static_assert(isSparsitySubsetLowerTriangle(SparsityStatic(Init22::sparsity),
                                              FactorData22::sparsity,
                                              FactorData22::permutation));
  static_assert(isSparsitySubset(
      SparsityStatic(Init21::sparsity), FactorData21::sparsity,
      FactorData21::permutation_row, FactorData21::permutation_col));
  constexpr auto num_rows = std::size_t{FactorData22::sparsity.numRows()};
  factorizeEntryWiseSubdiagonalImpl(factor21, init21, factor11,
                                    std::make_index_sequence<num_rows>());
  factorizeEntryWiseImpl(factor22, init22, factor21, factor11, regularization,
                         std::make_index_sequence<num_rows>());
}

/**
 * Factorize the given diagonal block \p factor using the initial values given
 * by \p init in an entry-wise way.
 *
 * Note that this implicitly assumes that no other blocks have influence on this
 * diagonal block, since no partial result can be returned that other
 * contributions could be added to. Instead the final result is computed in one
 * go directly starting with the input values.
 */
template <class FactorData, class Init>
void factorizeEntryWise(FactorData& factor, const Init& init,
                        const Regularization auto& regularization) {
  constexpr EmptyFactorDataLeft<FactorData> empty21;
  constexpr EmptyFactorDataDiagonal<decltype(empty21)> empty11;
  constexpr EmptyMatrixInput<FactorData::sparsity.numRows(), 0> empty_init21;
  factorizeEntryWise(empty11, empty_init21, init, empty21, factor,
                     regularization);
}

template <class FactorData11, class Init21, class Init22, class FactorData21,
          class FactorData22>
void factorize(const FactorData11& factor11, const Init21& init21,
               const Init22& init22, FactorData21& factor21,
               FactorData22& factor22,
               const Regularization auto& regularization,
               FactorizeMethodEntryWise) {
  factorizeEntryWise(factor11, init21, init22, factor21, factor22,
                     regularization);
}

template <class FactorData, class Init>
void factorize(FactorData& factor, const Init& init,
               const Regularization auto& regularization,
               FactorizeMethodEntryWise) {
  factorizeEntryWise(factor, init, regularization);
}

}  // namespace ctldl
