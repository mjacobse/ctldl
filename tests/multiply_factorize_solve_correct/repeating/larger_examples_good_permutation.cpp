#include "tests/multiply_factorize_solve_correct/repeating/test_functor.hpp"
#include "tests/test_matrices/repeating/nos4.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/test_set.hpp"
#include "tests/utility/test_set_foreach.hpp"

#include <boost/test/unit_test.hpp>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeating)
BOOST_AUTO_TEST_CASE(LargerExamplesGoodPermutation) {
  const auto matrix_permutation_pairs =
      (makeTypeArgument<TestMatrixNos4A<double>, TestMatrixNos4A<float>>() ^
       makeTypeArgument<TestMatrixNos4B<double>, TestMatrixNos4B<float>>()) *
      makeTypeArgument<TestPermutationNos4>();

  const auto factorize_value_types = makeTypeArgument<double, float>();
  const auto factorize_method =
      makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
  const auto solution_generators = makeValueArgument(
      {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
       getSolutionGeneratorNormallyDistributed(1000.0)});
  const auto repetition_counts = makeValueArgument({0, 1, 2, 9});

  const auto test_set = matrix_permutation_pairs * factorize_value_types *
                        factorize_method * solution_generators *
                        repetition_counts;
  foreach<TesterMultiplyFactorizeSolveCorrectRepeating>(test_set);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
