#include "tests/multiply_factorize_solve_correct/repeating/test_functor.hpp"
#include "tests/test_matrices/repeating/nos4.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/test_set.hpp"
#include "tests/utility/test_set_foreach.hpp"

#include <boost/test/unit_test.hpp>

#include <functional>
#include <random>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeating)
BOOST_AUTO_TEST_CASE(LargerExamplesGoodPermutation) {
  static constexpr auto matrix_permutation_pairs =
      (makeTypeArgument<TestMatrixNos4A<double>, TestMatrixNos4A<float>>() ^
       makeTypeArgument<TestMatrixNos4B<double>, TestMatrixNos4B<float>>()) *
      makeTypeArgument<TestPermutationNos4>();

  static constexpr auto factorize_value_types =
      makeTypeArgument<double, float>();
  static constexpr auto factorize_method =
      makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
  const auto solution_generators = makeValueArgument(
      {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
       getSolutionGeneratorNormallyDistributed(1000.0)});
  const auto repetition_counts = makeValueArgument<std::size_t>({0, 1, 2, 9});

  std::mt19937 value_generator{0};
  static constexpr auto test_set_types =
      matrix_permutation_pairs * factorize_value_types * factorize_method;
  const auto test_set_values = solution_generators * repetition_counts *
                               makeValueArgument({std::ref(value_generator)});

  foreach
    <TesterMultiplyFactorizeSolveCorrectRepeating, ^^test_set_types>(
        test_set_values);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
