#pragma once

#include <ctldl/sparsity/sparsity_csr.hpp>

#include <cstddef>
#include <vector>

namespace ctldl {

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
constexpr auto getInfluencedList(const std::size_t i, const std::size_t j,
                                 const SparsityViewCSR sparsity_source,
                                 const SparsityViewCSR sparsity_target,
                                 const std::size_t target_column_limit) {
  std::vector<Influenced> influenced_list;
  for (std::size_t k = 0; k < target_column_limit; ++k) {
    if (sparsity_source.isNonZero(k, j)) {
      influenced_list.push_back(Influenced{sparsity_source.entryIndex(k, j),
                                           sparsity_target.entryIndex(i, k)});
    }
  }
  return influenced_list;
}

constexpr auto getInfluencedList(const std::size_t i, const std::size_t j,
                                 const SparsityViewCSR sparsity_source,
                                 const SparsityViewCSR sparsity_target) {
  return getInfluencedList(i, j, sparsity_source, sparsity_target,
                           sparsity_target.numCols());
}

/**
 * Same as getInfluencedList(), but limits the influences to those that end up
 * in the lower triangle of the target sparsity.
 */
constexpr auto getInfluencedListLowerTriangle(
    const std::size_t i, const std::size_t j,
    const SparsityViewCSR sparsity_source,
    const SparsityViewCSR sparsity_target) {
  return getInfluencedList(i, j, sparsity_source, sparsity_target, i);
}

}  // namespace ctldl
