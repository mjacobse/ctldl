#pragma once

#include <ctldl/sparsity/sparsity.hpp>

#include <cstddef>
#include <vector>

namespace ctldl {

struct Contribution {
  std::size_t entry_index_ik;
  std::size_t entry_index_jk;
  std::size_t k;
};

template <class TernaryFunction>
constexpr void foreachContributingEntry(const SparsityViewCSR& sparsity,
                                        const SparsityViewCSR& sparsity_below,
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

constexpr auto getContributionsMixed(const SparsityViewCSR sparsity,
                                     const SparsityViewCSR sparsity_below,
                                     const std::size_t i, const std::size_t j,
                                     const std::size_t column_limit) {
  std::vector<Contribution> contributions;
  const auto add_contribution = [&contributions](
                                    const std::size_t entry_index_ik,
                                    const std::size_t entry_index_jk,
                                    const std::size_t col_index) {
    contributions.push_back({entry_index_ik, entry_index_jk, col_index});
  };
  foreachContributingEntry(sparsity, sparsity_below, i, j, column_limit,
                           add_contribution);
  return contributions;
}

constexpr auto getContributionsMixed(const SparsityViewCSR sparsity,
                                     const SparsityViewCSR sparsity_below,
                                     const std::size_t i, const std::size_t j) {
  return getContributionsMixed(
      sparsity, sparsity_below, i, j,
      std::min(sparsity.numCols(), sparsity_below.numCols()));
}

constexpr auto getContributionsLowerTriangle(const SparsityViewCSR sparsity,
                                             const std::size_t i,
                                             const std::size_t j) {
  return getContributionsMixed(sparsity, sparsity, i, j, j);
}

constexpr auto getContributionsRectangular(const SparsityViewCSR sparsity,
                                           const std::size_t i,
                                           const std::size_t j) {
  return getContributionsMixed(sparsity, sparsity, i, j);
}

}  // namespace ctldl
