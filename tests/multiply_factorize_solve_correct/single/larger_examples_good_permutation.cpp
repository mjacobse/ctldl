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
  static constexpr auto matrix_permutation_pairs =
      makeTypeArgument<TestMatrixLFAT5<double>, TestMatrixLFAT5<float>>() *
      makeTypeArgument<TestPermutationLFAT5>();

  static constexpr auto factorize_value_types =
      makeTypeArgument<double, float>();
  static constexpr auto factorize_method =
      makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
  const auto solution_generators = makeValueArgument(
      {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
       getSolutionGeneratorNormallyDistributed(1000.0)});

  std::mt19937 value_generator{0};
  static constexpr auto test_set_types =
      matrix_permutation_pairs * factorize_value_types * factorize_method;
  const auto test_set_values =
      solution_generators * makeValueArgument({std::ref(value_generator)});
  foreach
    <TesterMultiplyFactorizeSolveCorrectSingle, ^^test_set_types>(
        test_set_values);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
