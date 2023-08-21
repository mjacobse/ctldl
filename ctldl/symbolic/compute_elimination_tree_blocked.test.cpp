#include <ctldl/symbolic/compute_elimination_tree_blocked.hpp>

#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>

#include "tests/utility/test_static.hpp"

#include <boost/test/unit_test.hpp>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_CASE(ComputeEliminationTreeBlockedNos2Test) {
  constexpr auto sparsity11 = makeEmptySparsityCSR<3, 3>();
  constexpr auto sparsity21 =
      makeSparsityCSR<3, 3>({{0, 0}, {1, 1}, {1, 2}, {2, 1}, {2, 2}});
  constexpr auto sparsity22 = makeEmptySparsityCSR<3, 3>();
  constexpr auto tree =
      computeEliminationTreeBlocked(sparsity11, sparsity21, sparsity22);

  constexpr auto tree_expected = [] {
    EliminationTree<6> tree;
    tree.parent[0] = 3;
    tree.parent[1] = 4;
    tree.parent[2] = 4;
    tree.parent[4] = 5;
    return tree;
  }();

  CTLDL_TEST_STATIC(tree.parent == tree_expected.parent,
                    boost::test_tools::per_element());
}

}  // anonymous namespace
}  // namespace ctldl
