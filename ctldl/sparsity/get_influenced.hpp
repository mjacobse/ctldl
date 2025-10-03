#pragma once

#include <ctldl/utility/fix_init_if_zero_length_array.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <class SparsitySource>
constexpr std::size_t getNumInfluenced(const SparsitySource& sparsity_source,
                                       const std::size_t j,
                                       const std::size_t k_end) {
  std::size_t num = 0;
  for (std::size_t k = 0; k < k_end; ++k) {
    if (sparsity_source.isNonZero(k, j)) {
      num += 1;
    }
  }
  return num;
};

struct Influenced {
  std::size_t entry_index_source;
  std::size_t entry_index_target;
};

/**
 * Returns information of the influence that entries in block \p sparsity_source
 * have on the block \p sparsity_target together with the entry (i, j) in the
 * linking block.
 *
 * With S: \p sparsity_source, T: \p sparsity_target and L being the linking
 * block between the two, the situation looks like this:
 *
 * [*            ]
 * [* *          ]
 * [* * *        ]
 * [* S * *      ]
 * [* * * * *    ]
 * [* L * T * *  ]
 * [* * * * * * *]
 *
 * So this function will return a list of mappings from entries in the
 * j-th column of \p sparsity_source to corresponding entries in the i-th row of
 * \p sparsity_target that get influenced by that source entry together with the
 * entry (i,j) in the linking block.
 *
 * Only considers entries in the target sparsity that are in a column earlier
 * than \p target_column_limit.
 */
template <std::size_t i, std::size_t j, auto sparsity_source,
          auto sparsity_target,
          std::size_t target_column_limit = sparsity_target.numCols()>
constexpr auto getInfluencedList() {
  constexpr auto num_influenced =
      getNumInfluenced(sparsity_source, j, target_column_limit);
  std::array<Influenced, num_influenced> influenced_list;
  fixInitIfZeroLengthArray(influenced_list);
  std::size_t influenced_index = 0;
  for (std::size_t k = 0; k < target_column_limit; ++k) {
    if (sparsity_source.isNonZero(k, j)) {
      influenced_list[influenced_index] = Influenced{
          sparsity_source.entryIndex(k, j), sparsity_target.entryIndex(i, k)};
      influenced_index += 1;
    }
  }
  return influenced_list;
}

/**
 * Same as getInfluencedList(), but limits the influences to those that end up
 * in the lower triangle of the target sparsity.
 */
template <std::size_t i, std::size_t j, auto sparsity_source,
          auto sparsity_target>
constexpr auto getInfluencedListLowerTriangle() {
  return getInfluencedList<i, j, sparsity_source, sparsity_target, i>();
}

}  // namespace ctldl
