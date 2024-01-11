#pragma once

#include <ctldl/factor_data/empty_factor_data_diagonal.hpp>
#include <ctldl/factor_data/empty_factor_data_left.hpp>
#include <ctldl/factorize/factor_initialization.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/factorize/fill_with_original_matrix.hpp>
#include <ctldl/matrix/empty_matrix_input.hpp>
#include <ctldl/sparsity/get_influenced.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/sparsity/is_sparsity_subset.hpp>
#include <ctldl/symbolic/is_chordal_blocked.hpp>
#include <ctldl/utility/make_index_sequence.hpp>

#include <cstddef>
#include <utility>

namespace ctldl {

/**
 * Finishes factorization of the entry given by \p entry_index in the given
 * diagonal block \p fact .
 *
 * \return The contribution of the finished entry to the diagonal entry in its
 * row.
 */
template <std::size_t entry_index, class FactorData>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerSelf(
    FactorData& fact) {
  constexpr auto i = FactorData::sparsity.entries[entry_index].row_index;
  constexpr auto j = FactorData::sparsity.entries[entry_index].col_index;

  const auto Lij_scaled = fact.L[entry_index];

  static constexpr auto influenced_list =
      getInfluencedListLowerTriangle<i, j, FactorData::sparsity,
                                     FactorData::sparsity>();
  for (const auto influenced : influenced_list) {
    fact.L[influenced.entry_index_target] -=
        fact.L[influenced.entry_index_source] * Lij_scaled;
  }

  const auto Lij = Lij_scaled / fact.D[j];
  fact.L[entry_index] = Lij;
  return Lij_scaled * Lij;
}

template <std::size_t... EntryIndices, class FactorData>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerSelf(
    FactorData& fact, typename FactorData::Value Di_init,
    std::index_sequence<EntryIndices...>) {
  return (Di_init - ... - factorizeUpLookingInnerSelf<EntryIndices>(fact));
}

/**
 * Finishes factorization of the entry given by \p entry_index in block
 * \p factor21 (L21) after applying its contributions to entries to the right
 * in block \p factor22 (L22).
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
template <std::size_t entry_index, class FactorData11, class FactorData21,
          class FactorData22>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorData11& factor11, FactorData21& factor21,
    FactorData22& factor22) {
  constexpr auto i = FactorData21::sparsity.entries[entry_index].row_index;
  constexpr auto j = FactorData21::sparsity.entries[entry_index].col_index;

  const auto Lij_scaled = factor21.L[entry_index];

  static constexpr auto influenced_list21 =
      getInfluencedList<i, j, FactorData11::sparsity, FactorData21::sparsity>();
  for (const auto influenced : influenced_list21) {
    factor21.L[influenced.entry_index_target] -=
        factor11.L[influenced.entry_index_source] * Lij_scaled;
  }
  static constexpr auto influenced_list22 =
      getInfluencedListLowerTriangle<i, j, FactorData21::sparsity,
                                     FactorData22::sparsity>();
  for (const auto influenced : influenced_list22) {
    factor22.L[influenced.entry_index_target] -=
        factor21.L[influenced.entry_index_source] * Lij_scaled;
  }

  const auto Lij = Lij_scaled / factor11.D[j];
  factor21.L[entry_index] = Lij;
  return Lij_scaled * Lij;
}

/**
 * Finishes factorization of the entry given by \p entry_index in \p factor31
 * (L31) after applying its contributions to entries to the right in
 * \p factor31 (L31), \p factor32 (L32) and \p factor33 (L33).
 *
 * [*                  ]
 * [* L11              ]
 * [*  *  *            ]
 * [* L21 *  *         ]
 * [*  *  *  *  *      ]
 * [* L31 * L32 * L33  ]
 * [*  *  *  *  *  *  *]
 */
template <std::size_t entry_index, class FactorData11, class FactorData21,
          class FactorData31, class FactorData32, class FactorData33>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorData11& factor11, const FactorData21& factor21,
    FactorData31& factor31, FactorData32& factor32, FactorData33& factor33) {
  constexpr auto& sparsity11 = FactorData11::sparsity;
  constexpr auto& sparsity21 = FactorData21::sparsity;
  constexpr auto& sparsity31 = FactorData31::sparsity;
  constexpr auto& sparsity32 = FactorData32::sparsity;
  constexpr auto& sparsity33 = FactorData33::sparsity;

  constexpr auto i = sparsity31.entries[entry_index].row_index;
  constexpr auto j = sparsity31.entries[entry_index].col_index;

  const auto Lij_scaled = factor31.L[entry_index];

  static constexpr auto influenced_list31 =
      getInfluencedList<i, j, sparsity11, sparsity31>();
  for (const auto influenced : influenced_list31) {
    factor31.L[influenced.entry_index_target] -=
        factor11.L[influenced.entry_index_source] * Lij_scaled;
  }

  static constexpr auto influenced_list32 =
      getInfluencedList<i, j, sparsity21, sparsity32>();
  for (const auto influenced : influenced_list32) {
    factor32.L[influenced.entry_index_target] -=
        factor21.L[influenced.entry_index_source] * Lij_scaled;
  }

  static constexpr auto influenced_list33 =
      getInfluencedListLowerTriangle<i, j, sparsity31, sparsity33>();
  for (const auto influenced : influenced_list33) {
    factor33.L[influenced.entry_index_target] -=
        factor31.L[influenced.entry_index_source] * Lij_scaled;
  }

  const auto Lij = Lij_scaled / factor11.D[j];
  factor31.L[entry_index] = Lij;
  return Lij_scaled * Lij;
}

template <std::size_t... EntryIndices, class FactorData11, class FactorData21,
          class FactorData31, class FactorData32, class FactorData33>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorData11& factor11, const FactorData21& factor21,
    FactorData31& factor31, FactorData32& factor32, FactorData33& factor33,
    typename FactorData33::Value Di_init,
    std::index_sequence<EntryIndices...>) {
  return (Di_init - ... -
          factorizeUpLookingInnerLeft<EntryIndices>(
              factor11, factor21, factor31, factor32, factor33));
}

template <std::size_t... EntryIndices, class FactorData11, class FactorData21,
          class FactorData22>
[[gnu::always_inline]] inline auto factorizeUpLookingInnerLeft(
    const FactorData11& factor11, FactorData21& factor21,
    FactorData22& factor22, typename FactorData22::Value Di_init,
    std::index_sequence<EntryIndices...>) {
  return (
      Di_init - ... -
      factorizeUpLookingInnerLeft<EntryIndices>(factor11, factor21, factor22));
}

template <std::size_t i, class FactorData11, class Init21, class Init22,
          class FactorData21, class FactorData22>
[[gnu::always_inline]] inline void factorizeUpLookingImpl(
    const FactorData11& factor11, const Init21& init21, const Init22& init22,
    FactorData21& factor21, FactorData22& factor22) {
  constexpr auto& sparsity21 = FactorData21::sparsity;
  constexpr auto& sparsity22 = FactorData22::sparsity;
  static_assert(sparsity22.num_rows == sparsity21.num_rows);

  auto Di = getInitialFactorValueD<i>(init22, factor22);

  initializeFactorRow<i>(init21, factor21);
  initializeFactorRow<i>(init22, factor22);

  constexpr auto row_begin21 = sparsity21.row_begin_indices[i];
  constexpr auto row_end21 = sparsity21.row_begin_indices[i + 1];
  Di = factorizeUpLookingInnerLeft(factor11, factor21, factor22, Di,
                                   makeIndexSequence<row_begin21, row_end21>());
  constexpr auto row_begin22 = sparsity22.row_begin_indices[i];
  constexpr auto row_end22 = sparsity22.row_begin_indices[i + 1];
  Di = factorizeUpLookingInnerSelf(factor22, Di,
                                   makeIndexSequence<row_begin22, row_end22>());
  factor22.D[i] = Di;
}

template <std::size_t... RowIndices, class FactorData11, class Init21,
          class Init22, class FactorData21, class FactorData22>
void factorizeUpLookingImpl(const FactorData11& factor11, const Init21& init21,
                            const Init22& init22, FactorData21& factor21,
                            FactorData22& factor22,
                            std::index_sequence<RowIndices...>) {
  (factorizeUpLookingImpl<RowIndices>(factor11, init21, init22, factor21,
                                      factor22),
   ...);
}

/**
 * Factorize blocks \p factor21 (L21) and \p factor22 (L22) using the initial
 * values given by \p init21 and \p init22 for those blocks and the already
 * finished block \p factor22 (L11) in an up-looking way.
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
 * starting with the initial values.
 */
template <class FactorData11, class Init21, class Init22, class FactorData21,
          class FactorData22>
void factorizeUpLooking(const FactorData11& factor11, const Init21& init21,
                        const Init22& init22, FactorData21& factor21,
                        FactorData22& factor22) {
  static_assert(isChordalBlocked(FactorData11::sparsity, FactorData21::sparsity,
                                 FactorData22::sparsity));
  static_assert(isSparsitySubsetLowerTriangle<Init22::sparsity>(
      FactorData22::sparsity, FactorData22::permutation));
  static_assert(isSparsitySubset(Init21::sparsity, FactorData21::sparsity,
                                 FactorData21::permutation_row,
                                 FactorData21::permutation_col));
  constexpr auto num_rows = std::size_t{FactorData22::sparsity.num_rows};
  factorizeUpLookingImpl(factor11, init21, init22, factor21, factor22,
                         std::make_index_sequence<num_rows>());
}

/**
 * Factorize the given diagonal block \p factor using the initial values given
 * by \p init in an up-looking way.
 *
 * Note that this implicitly assumes that no other blocks have influence on this
 * diagonal block, since no partial result can be returned that other
 * contributions could be added to. Instead the final result is computed in one
 * go directly starting with the initial values.
 */
template <class FactorData, class Init>
void factorizeUpLooking(FactorData& factor, const Init& init) {
  constexpr EmptyFactorDataLeft<FactorData> empty21;
  constexpr EmptyFactorDataDiagonal<decltype(empty21)> empty11;
  constexpr EmptyMatrixInput<FactorData::sparsity.num_rows, 0> empty_init21;
  factorizeUpLooking(empty11, empty_init21, init, empty21, factor);
}

template <class FactorData11, class Init21, class Init22, class FactorData21,
          class FactorData22>
void factorize(const FactorData11& factor11, const Init21& init21,
               const Init22& init22, FactorData21& factor21,
               FactorData22& factor22, FactorizeMethodUpLooking) {
  factorizeUpLooking(factor11, init21, init22, factor21, factor22);
}

template <class FactorData, class Init>
void factorize(FactorData& factor, const Init& init, FactorizeMethodUpLooking) {
  factorizeUpLooking(factor, init);
}

template <std::size_t i, class FactorData11, class FactorData21, class Init31,
          class Init32, class Init33, class FactorData31, class FactorData32,
          class FactorData33>
void factorizePartialUpLookingImpl(const FactorData11& factor11,
                                   const FactorData21& factor21,
                                   const Init31& init31, const Init32& init32,
                                   const Init33& init33, FactorData31& factor31,
                                   FactorData32& factor32,
                                   FactorData33& factor33) {
  constexpr auto& sparsity31 = FactorData31::sparsity;

  initializeFactorRow<i>(init31, factor31);
  initializeFactorRow<i>(init32, factor32);
  initializeFactorRow<i>(init33, factor33);
  auto Di = factor33.D[i];
  constexpr auto row_begin31 = sparsity31.row_begin_indices[i];
  constexpr auto row_end31 = sparsity31.row_begin_indices[i + 1];
  Di = factorizeUpLookingInnerLeft(factor11, factor21, factor31, factor32,
                                   factor33, Di,
                                   makeIndexSequence<row_begin31, row_end31>());
  factor33.D[i] = Di;
}

template <std::size_t... RowIndices, class FactorData11, class FactorData21,
          class Init31, class Init32, class Init33, class FactorData31,
          class FactorData32, class FactorData33>
void factorizePartialUpLookingImpl(const FactorData11& factor11,
                                   const FactorData21& factor21,
                                   const Init31& init31, const Init32& init32,
                                   const Init33& init33, FactorData31& factor31,
                                   FactorData32& factor32,
                                   FactorData33& factor33,
                                   std::index_sequence<RowIndices...>) {
  (factorizePartialUpLookingImpl<RowIndices>(factor11, factor21, init31, init32,
                                             init33, factor31, factor32,
                                             factor33),
   ...);
}

/**
 * Finish factorization of \p factor31 (L31) using already finished blocks
 * \p factor21 (L21) and \p factor11 (L11). Also apply influences from
 * \p factor31 (L31) onto \p factor32 (L32) and \p factor33 (L33). The blocks
 * L32 and L33 are then not factorized further though, yielding only a partial
 * factorization for them.
 *
 * [*                  ]
 * [* L11              ]
 * [*  *  *            ]
 * [* L21 *  *         ]
 * [*  *  *  *  *      ]
 * [* L31 * L32 * L33  ]
 * [*  *  *  *  *  *  *]
 *
 * L31, L32 and L33 are initialized with the values given by \p init31 ,
 * \p init32 and \p init33 .
 */
template <class FactorData11, class FactorData21, class Init31, class Init32,
          class Init33, class FactorData31, class FactorData32,
          class FactorData33>
void factorizePartialUpLooking(const FactorData11& factor11,
                               const FactorData21& factor21,
                               const Init31& init31, const Init32& init32,
                               const Init33& init33, FactorData31& factor31,
                               FactorData32& factor32, FactorData33& factor33) {
  constexpr auto num_rows = std::size_t{FactorData31::sparsity.num_rows};
  factorizePartialUpLookingImpl(factor11, factor21, init31, init32, init33,
                                factor31, factor32, factor33,
                                std::make_index_sequence<num_rows>());
}

}  // namespace ctldl
