#pragma once

#include <algorithm>
#include <array>
#include <cstddef>

namespace ctldl {

struct Contribution {
  std::size_t entry_index_ik;
  std::size_t entry_index_jk;
  std::size_t k;
};

template <class Sparsity, class SparsityBelow>
constexpr auto getNumContributionsMixed(const std::size_t i,
                                        const std::size_t j,
                                        const std::size_t column_limit) {
  std::size_t num = 0;
  for (std::size_t k = 0; k < column_limit; ++k) {
    num += SparsityBelow::isNonZero(i, k) && Sparsity::isNonZero(j, k);
  }
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
  for (std::size_t k = 0; k < column_limit; ++k) {
    if (SparsityBelow::isNonZero(i, k) && Sparsity::isNonZero(j, k)) {
      const std::size_t entry_index_ik = SparsityBelow::entryIndex(i, k);
      const std::size_t entry_index_jk = Sparsity::entryIndex(j, k);
      contributions[contrib_index] = {entry_index_ik, entry_index_jk, k};
      contrib_index += 1;
    }
  }
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
