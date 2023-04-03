#pragma once

#include <algorithm>
#include <array>

namespace ctldl {

struct Contribution {
  int entry_index_ik;
  int entry_index_jk;
  int k;
};

template <class Sparsity, class SparsityBelow>
constexpr auto getNumContributionsMixed(const int i, const int j,
                                        const int column_limit) {
  int num = 0;
  for (int k = 0; k < column_limit; ++k) {
    num += SparsityBelow::is_nonzero[i][k] && Sparsity::is_nonzero[j][k];
  }
  return num;
};

template <class Sparsity, class SparsityBelow, int i, int j,
          int column_limit = std::min(Sparsity::num_cols,
                                      SparsityBelow::num_cols)>
constexpr auto getContributionsMixed() {
  constexpr int num_contrib =
      getNumContributionsMixed<Sparsity, SparsityBelow>(i, j, column_limit);
  std::array<Contribution, num_contrib> contributions{};
  int contrib_index = 0;
  for (int k = 0; k < column_limit; ++k) {
    if (SparsityBelow::is_nonzero[i][k] && Sparsity::is_nonzero[j][k]) {
      const auto entry_index_ik = SparsityBelow::entryIndex(i, k);
      const auto entry_index_jk = Sparsity::entryIndex(j, k);
      contributions[contrib_index] = {entry_index_ik, entry_index_jk, k};
      contrib_index += 1;
    }
  }
  return contributions;
};

template <class Sparsity, int i, int j>
constexpr auto getContributionsLowerTriangle() {
  return getContributionsMixed<Sparsity, Sparsity, i, j, j>();
};

template <class Sparsity, int i, int j>
constexpr auto getContributionsRectangular() {
  return getContributionsMixed<Sparsity, Sparsity, i, j>();
};

}  // namespace ctldl
