#pragma once

#include <ctldl/empty_factor_data_left.hpp>
#include <ctldl/permutation/permuted_entry_lower_triangle.hpp>
#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/sparsity/get_matrix_value_at.hpp>
#include <ctldl/utility/make_index_sequence.hpp>
#include <ctldl/utility/square.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>

namespace ctldl {

template <std::size_t i, class FactorData>
[[gnu::always_inline]] inline auto applyContributionsRowDiagonal(
    const FactorData& left, const typename FactorData::Value value_init) {
  using Sparsity = typename FactorData::Sparsity;

  constexpr auto row_begin = std::size_t{Sparsity::row_begin_indices[i]};
  constexpr auto row_end = std::size_t{Sparsity::row_begin_indices[i + 1]};
  auto value = value_init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const std::size_t j = Sparsity::entries[entry_index_ij].col_index;
    value -= square(left.L[entry_index_ij]) * left.D[j];
  }
  return value;
}

template <class FactorData, std::size_t num_contributions>
[[gnu::always_inline]] inline auto applyContributions(
    const FactorData& fact,
    const std::array<Contribution, num_contributions>& contributions,
    const typename FactorData::Value value_init) {
  auto value = value_init;
  for (const auto c : contributions) {
    value -= fact.L[c.entry_index_ik] * fact.L[c.entry_index_jk] * fact.D[c.k];
  }
  return value;
}

template <std::size_t entry_index_ij, class FactorData, class Matrix,
          class FactorDataLeft>
[[gnu::always_inline]] inline auto factorizeImplRow(
    FactorData& self, const Matrix& input, const FactorDataLeft& left) {
  using Sparsity = typename FactorData::Sparsity;
  using SparsityLeft = typename FactorDataLeft::Sparsity;

  constexpr auto i = std::size_t{Sparsity::entries[entry_index_ij].row_index};
  constexpr auto j = std::size_t{Sparsity::entries[entry_index_ij].col_index};

  constexpr auto entry_orig =
      permutedEntryLowerTriangle(Entry{i, j}, FactorData::permutation);
  auto Lij =
      getMatrixValueAt<entry_orig.row_index, entry_orig.col_index>(input);

  static constexpr auto contributions_left =
      getContributionsRectangular<SparsityLeft, i, j>();
  static constexpr auto contributions_self =
      getContributionsLowerTriangle<Sparsity, i, j>();

  Lij = applyContributions(left, contributions_left, Lij);
  Lij = applyContributions(self, contributions_self, Lij);
  Lij /= self.D[j];
  self.L[entry_index_ij] = Lij;

  return square(Lij) * self.D[j];
}

template <std::size_t... EntryIndices, class FactorData, class Matrix,
          class FactorDataLeft>
[[gnu::always_inline]] inline auto factorizeImplRow(
    FactorData& self, const Matrix& input, const FactorDataLeft& left,
    const typename FactorData::Value Di_init,
    std::index_sequence<EntryIndices...>) {
  return (Di_init - ... - factorizeImplRow<EntryIndices>(self, input, left));
}

template <std::size_t i, class FactorData, class Matrix,
          class FactorDataLeft>
[[gnu::always_inline]] inline void factorizeImpl(FactorData& self,
                                                 const Matrix& input,
                                                 const FactorDataLeft& left) {
  using Sparsity = typename FactorData::Sparsity;
  using SparsityLeft = typename FactorDataLeft::Sparsity;
  static_assert(Sparsity::num_rows == SparsityLeft::num_rows);

  constexpr auto row_begin = std::size_t{Sparsity::row_begin_indices[i]};
  constexpr auto row_end = std::size_t{Sparsity::row_begin_indices[i + 1]};

  constexpr auto i_orig = FactorData::permutation[i];
  auto Di = getMatrixValueAt<i_orig, i_orig>(input);
  Di = applyContributionsRowDiagonal<i>(left, Di);
  Di = factorizeImplRow(self, input, left, Di,
                        makeIndexSequence<row_begin, row_end>());
  self.D[i] = Di;
}

template <std::size_t... RowIndices, class FactorData, class Matrix,
          class FactorDataLeft>
void factorizeImpl(FactorData& self, const Matrix& input,
                   const FactorDataLeft& left,
                   std::index_sequence<RowIndices...>) {
  (factorizeImpl<RowIndices>(self, input, left), ...);
}

template <class FactorData, class Matrix, class FactorDataLeft>
void factorize(FactorData& self, const Matrix& input,
               const FactorDataLeft& left) {
  constexpr auto num_rows = std::size_t{FactorData::Sparsity::num_rows};
  factorizeImpl(self, input, left, std::make_index_sequence<num_rows>());
}

template <class FactorData, class Matrix>
void factorize(FactorData& self, const Matrix& input) {
  constexpr EmptyFactorDataLeft<FactorData> empty_left;
  factorize(self, input, empty_left);
}

}  // namespace ctldl
