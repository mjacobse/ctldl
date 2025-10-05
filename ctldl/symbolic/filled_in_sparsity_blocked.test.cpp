#include <ctldl/symbolic/filled_in_sparsity_blocked.hpp>

#include <ctldl/sparsity/is_sparsity_equal.hpp>
#include <ctldl/sparsity/sparsity.hpp>

#include "tests/utility/test_sparsity_equal.hpp"
#include "tests/utility/test_static.hpp"

#include <boost/test/unit_test.hpp>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_CASE(FilledInSparsityBlockedEmptyTest) {
  constexpr auto empty = makeEmptySparsityStatic<1, 1>();
  constexpr auto filled = getFilledInSparsityBlocked<empty, empty, empty>();
  CTLDL_TEST_SPARSITY_EQUAL(filled.block11, empty);
  CTLDL_TEST_SPARSITY_EQUAL(filled.block21, empty);
  CTLDL_TEST_SPARSITY_EQUAL(filled.block22, empty);
}

BOOST_AUTO_TEST_CASE(FilledInSparsityBlockedNos2Test) {
  constexpr auto sparsity11 = makeEmptySparsityStatic<3, 3>();
  constexpr auto sparsity21 =
      makeSparsityStatic<3, 3>({{0, 0}, {1, 1}, {1, 2}, {2, 1}, {2, 2}});
  constexpr auto sparsity22 = sparsity11;

  constexpr auto filled =
      getFilledInSparsityBlocked<sparsity11, sparsity21, sparsity22>();

  constexpr auto expected = LowerTriangleBlocked{
      sparsity11, sparsity21, makeSparsityStatic<3, 3>({{2, 1}})};
  CTLDL_TEST_SPARSITY_EQUAL(filled.block11, expected.block11);
  CTLDL_TEST_SPARSITY_EQUAL(filled.block21, expected.block21);
  CTLDL_TEST_SPARSITY_EQUAL(filled.block22, expected.block22);
}

}  // anonymous namespace
}  // namespace ctldl
