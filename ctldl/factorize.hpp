#pragma once

#include <ctldl/sparsity/entry.hpp>
#include <ctldl/sparsity/get_contributions.hpp>
#include <ctldl/utility/square.hpp>

#include <array>
#include <cassert>

namespace ctldl {

template <int i, class FactorData>
[[gnu::always_inline]] inline auto applyContributionsRowDiagonal(
    const FactorData& left, const typename FactorData::Value value_init) {
  using Sparsity = typename FactorData::Sparsity;

  constexpr auto row_begin = Sparsity::row_begin_indices[i];
  constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
  auto value = value_init;
  for (auto entry_index_ij = row_begin; entry_index_ij != row_end;
       ++entry_index_ij) {
    const auto j = Sparsity::entries[entry_index_ij].col_index;
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

template <int i, int entry_index_ij, class FactorData, class Matrix,
          class FactorDataLeft>
[[gnu::always_inline]] inline auto factorizeImplRow(
    FactorData& self, const Matrix& input, const FactorDataLeft& left,
    const typename FactorData::Value Di_init) {
  using Sparsity = typename FactorData::Sparsity;
  using SparsityLeft = typename FactorDataLeft::Sparsity;

  constexpr auto row_end = Sparsity::row_begin_indices[i + 1];
  if constexpr (entry_index_ij < row_end) {
    constexpr auto j = Sparsity::entries[entry_index_ij].col_index;

    static constexpr auto contributions_left =
        getContributionsRectangular<SparsityLeft, i, j>();
    static constexpr auto contributions_self =
        getContributionsLowerTriangle<Sparsity, i, j>();

    auto Lij = input.get(i, j);
    Lij = applyContributions(left, contributions_left, Lij);
    Lij = applyContributions(self, contributions_self, Lij);
    Lij /= self.D[j];
    self.L[entry_index_ij] = Lij;

    const auto Di =
        factorizeImplRow<i, entry_index_ij + 1>(self, input, left, Di_init);
    return Di - square(Lij) * self.D[j];
  }
  return Di_init;
}

template <int i = 0, class FactorData, class Matrix, class FactorDataLeft>
[[gnu::always_inline]] inline void factorizeImpl(
    FactorData& self, const Matrix& input, const FactorDataLeft& left) {
  using Sparsity = typename FactorData::Sparsity;
  using SparsityLeft = typename FactorDataLeft::Sparsity;
  static_assert(Sparsity::num_rows == SparsityLeft::num_rows);

  if constexpr (i < Sparsity::num_rows) {
    auto Di = input.get(i, i);
    Di = applyContributionsRowDiagonal<i>(left, Di);
    Di = factorizeImplRow<i, Sparsity::row_begin_indices[i]>(
        self, input, left, Di);
    self.D[i] = Di;
    factorizeImpl<i + 1>(self, input, left);
  }
}

template <class FactorData, class Matrix, class FactorDataLeft>
void factorize(FactorData& self, const Matrix& input,
               const FactorDataLeft& left) {
  factorizeImpl(self, input, left);
}

template <int num_rows_>
struct EmptyLeftSparsity {
  static constexpr auto num_rows = num_rows_;
  static constexpr auto num_cols = 0;
  static constexpr std::array<int, num_rows + 1> row_begin_indices{};
  static constexpr std::array<Entry, 0> entries{};
  static constexpr std::array<std::array<bool, 0>, num_rows> is_nonzero{};

  static constexpr int entryIndex(const int /*i*/, const int /*j*/) {
    assert(false);
    return -1;
  }
};

template <int num_rows_, class Value_>
struct EmptyLeftFactorData {
  static constexpr auto num_rows = num_rows_;
  using Value = Value_;
  using Sparsity = EmptyLeftSparsity<num_rows>;
  static constexpr std::array<Value, 0> L{};
  static constexpr std::array<Value, 0> D{};
};

template <class FactorData, class Matrix>
void factorize(FactorData& self, const Matrix& input) {
  const EmptyLeftFactorData<FactorData::Sparsity::num_rows,
                            typename FactorData::Value>
      empty_left;
  factorize(self, input, empty_left);
}

}  // namespace ctldl
