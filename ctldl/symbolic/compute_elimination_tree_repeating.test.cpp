#include <ctldl/symbolic/compute_elimination_tree_repeating.hpp>

#include <ctldl/sparsity/sparsity_csr.hpp>

#include "tests/utility/test_static.hpp"

#include <boost/test/unit_test.hpp>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_CASE(ComputeEliminationTreeRepeatingExample0) {
  constexpr auto sparsity_A = makeSparsityCSR<5, 5>({{3, 2}});
  constexpr auto sparsity_B =
      makeSparsityCSR<5, 5>({{0, 3}, {2, 0}, {2, 1}, {2, 4}, {4, 4}});
  constexpr auto tree = computeEliminationTreeRepeating(sparsity_A, sparsity_B);

  constexpr auto tree_expected = [] {
    EliminationTree<2 * 5> ret;
    ret.parent[0] = 2;
    ret.parent[1] = 2 + 5;
    ret.parent[2] = 3;
    ret.parent[3] = 4;
    ret.parent[4] = 0 + 5;
    ret.parent[0 + 5] = 2 + 5;
    ret.parent[2 + 5] = 3 + 5;
    ret.parent[3 + 5] = 4 + 5;
    return ret;
  }();

  CTLDL_TEST_STATIC(tree.parent == tree_expected.parent,
                    boost::test_tools::per_element());
}

BOOST_AUTO_TEST_CASE(ComputeEliminationTreeRepeatingExample1) {
  constexpr auto sparsity_A = makeSparsityCSR<7, 7>({{1, 0}, {2, 0}, {2, 1}});
  constexpr auto sparsity_B = makeSparsityCSR<7, 7>(
      {{0, 5}, {0, 6}, {1, 3}, {1, 4}, {1, 5}, {2, 4}, {2, 6}});
  constexpr auto tree = computeEliminationTreeRepeating(sparsity_A, sparsity_B);

  constexpr auto tree_expected = [] {
    EliminationTree<2 * 7> ret;
    ret.parent[0] = 1;
    ret.parent[1] = 2;
    ret.parent[3] = 1 + 7;
    ret.parent[4] = 1 + 7;
    ret.parent[0 + 5] = 0 + 7;
    ret.parent[1 + 5] = 0 + 7;
    ret.parent[2 + 5] = 1 + 7;
    ret.parent[3 + 5] = 2 + 7;
    return ret;
  }();

  CTLDL_TEST_STATIC(tree.parent == tree_expected.parent,
                    boost::test_tools::per_element());
}

}  // anonymous namespace
}  // namespace ctldl
