#include <ctldl/symbolic/compute_elimination_tree_blocked.hpp>

#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_csr.hpp>

#include <boost/test/unit_test.hpp>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_CASE(ComputeEliminationTreeBlockedNos2Test) {
  constexpr auto sparsity11 = makeEmptySparsityStaticCSR<3, 3>();
  constexpr auto sparsity21 =
      makeSparsityStaticCSR<3, 3>({{0, 0}, {1, 1}, {1, 2}, {2, 1}, {2, 2}});
  constexpr auto sparsity22 = makeEmptySparsityStaticCSR<3, 3>();
  const auto tree =
      computeEliminationTreeBlocked(sparsity11, sparsity21, sparsity22);

  const auto tree_expected = [] {
    EliminationTree ret(6);
    ret.parent[0] = 3;
    ret.parent[1] = 4;
    ret.parent[2] = 4;
    ret.parent[4] = 5;
    return ret;
  }();

  BOOST_TEST(tree.parent == tree_expected.parent,
             boost::test_tools::per_element());
}

}  // anonymous namespace
}  // namespace ctldl
