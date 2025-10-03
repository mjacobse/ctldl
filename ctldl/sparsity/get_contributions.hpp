#pragma once

#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

struct Contribution {
  std::size_t entry_index_ik;
  std::size_t entry_index_jk;
  std::size_t k;
};

template <class Sparsity, class SparsityBelow, class TernaryFunction>
constexpr void foreachContributingEntry(const Sparsity& sparsity,
                                        const SparsityBelow& sparsity_below,
                                        const std::size_t i,
                                        const std::size_t j,
                                        const std::size_t column_limit,
                                        TernaryFunction f) {
  const auto row_begin_i = sparsity_below.rowBeginIndices()[i];
  const auto row_end_i = sparsity_below.rowBeginIndices()[i + 1];
  const auto row_begin_j = sparsity.rowBeginIndices()[j];
  const auto row_end_j = sparsity.rowBeginIndices()[j + 1];

  auto entry_index_ik = row_begin_i;
  auto entry_index_jk = row_begin_j;
  while (entry_index_ik != row_end_i && entry_index_jk != row_end_j) {
    const auto col_index_i = sparsity_below.entries()[entry_index_ik].col_index;
    const auto col_index_j = sparsity.entries()[entry_index_jk].col_index;
    if (col_index_i >= column_limit || col_index_j >= column_limit) {
      break;
    }
    if (col_index_i == col_index_j) {
      f(entry_index_ik, entry_index_jk, col_index_i);
    }
    entry_index_ik += (col_index_i <= col_index_j);
    entry_index_jk += (col_index_j <= col_index_i);
  }
}

template <class Sparsity, class SparsityBelow>
constexpr auto getNumContributionsMixed(const Sparsity& sparsity,
                                        const SparsityBelow& sparsity_below,
                                        const std::size_t i,
                                        const std::size_t j,
                                        const std::size_t column_limit) {
  std::size_t num = 0;
  const auto count_contribution = [&num](const std::size_t /*entry_index_ik*/,
                                         const std::size_t /*entry_index_jk*/,
                                         const std::size_t /*col_index*/) {
    num += 1;
  };
  foreachContributingEntry(sparsity, sparsity_below, i, j, column_limit,
                           count_contribution);
  return num;
}

template <auto sparsity, auto sparsity_below, std::size_t i, std::size_t j,
          std::size_t column_limit = std::min(sparsity.numCols(),
                                              sparsity_below.numCols())>
constexpr auto getContributionsMixed() {
  constexpr auto num_contrib =
      getNumContributionsMixed(sparsity, sparsity_below, i, j, column_limit);
  std::array<Contribution, num_contrib> contributions;
  fixInitIfZeroLengthArray(contributions);
  std::size_t contrib_index = 0;
  const auto add_contribution = [&contributions, &contrib_index](
                                    const std::size_t entry_index_ik,
                                    const std::size_t entry_index_jk,
                                    const std::size_t col_index) {
    contributions[contrib_index] = {entry_index_ik, entry_index_jk, col_index};
    contrib_index += 1;
  };
  foreachContributingEntry(sparsity, sparsity_below, i, j, column_limit,
                           add_contribution);
  return contributions;
}

template <auto sparsity, std::size_t i, std::size_t j>
constexpr auto getContributionsLowerTriangle() {
  return getContributionsMixed<sparsity, sparsity, i, j, j>();
}

template <auto sparsity, std::size_t i, std::size_t j>
constexpr auto getContributionsRectangular() {
  return getContributionsMixed<sparsity, sparsity, i, j>();
}

}  // namespace ctldl
