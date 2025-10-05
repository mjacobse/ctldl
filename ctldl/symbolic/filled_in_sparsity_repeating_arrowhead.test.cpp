#include <ctldl/symbolic/filled_in_sparsity_repeating_arrowhead.hpp>

#include <ctldl/sparsity/is_sparsity_equal.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include "tests/utility/test_sparsity_equal.hpp"
#include "tests/utility/test_static.hpp"

#include <boost/test/unit_test.hpp>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_CASE(FilledInSparsityRepeatingArrowheadEmptyTest) {
  constexpr auto empty = makeEmptySparsityStatic<1, 1>();
  constexpr PermutationStatic<1> permutation{};
  constexpr auto filled =
      getFilledInSparsityRepeatingArrowhead<empty, empty, empty, permutation,
                                            permutation>();
  CTLDL_TEST_SPARSITY_EQUAL(filled.outer, empty);
}

BOOST_AUTO_TEST_CASE(FilledInSparsityRepeatingArrowheadFromDiagonalTest) {
  constexpr auto sparsity_A = makeSparsityStatic<3, 3>({{1, 0}, {2, 1}});
  constexpr auto sparsity_B = makeEmptySparsityStatic<3, 3>();
  constexpr auto sparsity_C = makeSparsityStatic<1, 3>({{0, 0}});
  constexpr PermutationStatic<3> permutation{};
  constexpr PermutationStatic<1> permutation_outer{};
  constexpr auto filled =
      getFilledInSparsityRepeatingArrowhead<sparsity_A, sparsity_B, sparsity_C,
                                            permutation, permutation_outer>();

  constexpr auto filled_correct =
      makeSparsityStatic<1, 3>({{0, 0}, {0, 1}, {0, 2}});
  CTLDL_TEST_SPARSITY_EQUAL(filled.outer, filled_correct);
}

BOOST_AUTO_TEST_CASE(FilledInSparsityRepeatingArrowheadFromBothTest) {
  constexpr auto sparsity_A = makeSparsityStatic<3, 3>({{2, 1}});
  constexpr auto sparsity_B = makeSparsityStatic<3, 3>({{1, 0}});
  constexpr auto sparsity_C = makeSparsityStatic<1, 3>({{0, 0}});
  constexpr PermutationStatic<3> permutation{};
  constexpr PermutationStatic<1> permutation_outer{};
  constexpr auto filled =
      getFilledInSparsityRepeatingArrowhead<sparsity_A, sparsity_B, sparsity_C,
                                            permutation, permutation_outer>();

  constexpr auto filled_correct =
      makeSparsityStatic<1, 3>({{0, 0}, {0, 1}, {0, 2}});
  CTLDL_TEST_SPARSITY_EQUAL(filled.outer, filled_correct);
}

BOOST_AUTO_TEST_CASE(FilledInSparsityRepeatingArrowheadBackwardsTest) {
  constexpr auto sparsity_A = makeEmptySparsityStatic<4, 4>();
  constexpr auto sparsity_B =
      makeSparsityStatic<4, 4>({{0, 1}, {1, 2}, {2, 3}});
  constexpr auto sparsity_C = makeSparsityStatic<1, 4>({{0, 3}});
  constexpr PermutationStatic<4> permutation{};
  constexpr PermutationStatic<1> permutation_outer{};
  constexpr auto filled =
      getFilledInSparsityRepeatingArrowhead<sparsity_A, sparsity_B, sparsity_C,
                                            permutation, permutation_outer>();

  constexpr auto filled_correct =
      makeSparsityStatic<1, 4>({{0, 0}, {0, 1}, {0, 2}, {0, 3}});
  CTLDL_TEST_SPARSITY_EQUAL(filled.outer, filled_correct);
}

}  // anonymous namespace
}  // namespace ctldl
