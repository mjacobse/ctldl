#pragma once

#include <array>
#include <cstddef>

namespace ctldl {

struct Contribution {
  std::size_t entry_index_ik;
  std::size_t entry_index_jk;
  std::size_t k;
};

template <class Sparsity, class SparsityBelow, class TernaryFunction>
constexpr void foreachContributingEntry(const std::size_t i,
                                        const std::size_t j,
                                        const std::size_t column_limit,
                                        TernaryFunction f) {
  const auto row_begin_i = SparsityBelow::row_begin_indices[i];
  const auto row_end_i = SparsityBelow::row_begin_indices[i + 1];
  const auto row_begin_j = Sparsity::row_begin_indices[j];
  const auto row_end_j = Sparsity::row_begin_indices[j + 1];

  auto entry_index_ik = row_begin_i;
  auto entry_index_jk = row_begin_j;
  while (entry_index_ik != row_end_i && entry_index_jk != row_end_j) {
    const auto col_index_i = SparsityBelow::entries[entry_index_ik].col_index;
    const auto col_index_j = Sparsity::entries[entry_index_jk].col_index;
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
constexpr auto getNumContributionsMixed(const std::size_t i,
                                        const std::size_t j,
                                        const std::size_t column_limit) {
  std::size_t num = 0;
  const auto count_contribution = [&num](const std::size_t /*entry_index_ik*/,
                                         const std::size_t /*entry_index_jk*/,
                                         const std::size_t /*col_index*/) {
    num += 1;
  };
  foreachContributingEntry<Sparsity, SparsityBelow>(i, j, column_limit,
                                                    count_contribution);
  return num;
};

template <class Sparsity, class SparsityBelow, std::size_t i, std::size_t j,
          std::size_t column_limit = std::min(Sparsity::num_cols,
                                              SparsityBelow::num_cols)>
constexpr auto getContributionsMixed() {
  constexpr auto num_contrib =
      getNumContributionsMixed<Sparsity, SparsityBelow>(i, j, column_limit);
  std::array<Contribution, num_contrib> contributions{};
  std::size_t contrib_index = 0;
  const auto add_contribution = [&contributions, &contrib_index](
                                    const std::size_t entry_index_ik,
                                    const std::size_t entry_index_jk,
                                    const std::size_t col_index) {
    contributions[contrib_index] = {entry_index_ik, entry_index_jk, col_index};
    contrib_index += 1;
  };
  foreachContributingEntry<Sparsity, SparsityBelow>(i, j, column_limit,
                                                    add_contribution);
  return contributions;
};

template <class Sparsity, std::size_t i, std::size_t j>
constexpr auto getContributionsLowerTriangle() {
  return getContributionsMixed<Sparsity, Sparsity, i, j, j>();
};

template <class Sparsity, std::size_t i, std::size_t j>
constexpr auto getContributionsRectangular() {
  return getContributionsMixed<Sparsity, Sparsity, i, j>();
};

}  // namespace ctldl
