#include "tests/multiply_factorize_solve_correct/single/test_functor.hpp"
#include "tests/test_matrices/single/lfat5.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/test_set.hpp"
#include "tests/utility/test_set_foreach.hpp"

#include <ctldl/factorize/factorize_method.hpp>

#include <boost/test/unit_test.hpp>

#include <functional>
#include <random>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectSingle)

BOOST_AUTO_TEST_CASE(LargerExamplesGoodPermutation) {
  const auto matrix_permutation_pairs =
      makeTypeArgument<TestMatrixLFAT5<double>, TestMatrixLFAT5<float>>() *
      makeTypeArgument<TestPermutationLFAT5>();

  const auto factorize_value_types = makeTypeArgument<double, float>();
  const auto factorize_method =
      makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
  const auto solution_generators = makeValueArgument(
      {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
       getSolutionGeneratorNormallyDistributed(1000.0)});

  std::mt19937 value_generator{0};
  const auto test_set = matrix_permutation_pairs * factorize_value_types *
                        factorize_method * solution_generators *
                        makeValueArgument({std::ref(value_generator)});
  foreach<TesterMultiplyFactorizeSolveCorrectSingle>(test_set);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
