#pragma once

#include <ctldl/sparsity/entry_print.hpp>
#include <ctldl/sparsity/sort_entries_row_major_sorted.hpp>

#include "tests/utility/test_static.hpp"

#include <boost/test/unit_test.hpp>

#define CTLDL_TEST_SPARSITY_EQUAL(sparsity_lhs, sparsity_rhs)         \
  CTLDL_TEST_STATIC(isSparsityEqual((sparsity_lhs), (sparsity_rhs))); \
  {                                                                   \
    const auto sortedEntries = [](auto entries) {                     \
      sortEntriesRowMajorSorted(entries);                             \
      return entries;                                                 \
    };                                                                \
    CTLDL_TEST_STATIC(sortedEntries((sparsity_lhs).entries()) ==      \
                          sortedEntries((sparsity_rhs).entries()),    \
                      boost::test_tools::per_element());              \
  }
