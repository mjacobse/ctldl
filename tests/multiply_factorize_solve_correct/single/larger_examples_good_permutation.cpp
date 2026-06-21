#include "tests/multiply_factorize_solve_correct/single/test_functor.hpp"
#include "tests/test_matrices/single/lfat5.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/test_set.hpp"
#include "tests/utility/test_set_foreach.hpp"

#include <ctldl/factorize/factorize_method.hpp>

#include <boost/test/unit_test.hpp>

#include <functional>
#include <random>
#include <vector>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectSingle)

BOOST_AUTO_TEST_CASE(LargerExamplesGoodPermutation) {
  static constexpr auto matrix_permutation_pairs =
      makeTemplateArgument<TestMatrixLFAT5<double>, TestMatrixLFAT5<float>>() *
      makeTemplateArgument(std::vector{test_permutation_lfat5});

  static constexpr auto factorize_value_types =
      makeTemplateArgument<double, float>();
  static constexpr auto factorize_method =
      makeTemplateArgument<FactorizeMethodUpLooking,
                           FactorizeMethodEntryWise>();
  const auto solution_generators = makeFunctionArgument(
      {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
       getSolutionGeneratorNormallyDistributed(1000.0)});

  std::mt19937 value_generator{0};
  static constexpr auto test_set_template =
      matrix_permutation_pairs * factorize_value_types * factorize_method;
  const auto test_set_function =
      solution_generators * makeFunctionArgument({std::ref(value_generator)});
  foreach
    <^^TesterMultiplyFactorizeSolveCorrectSingle, ^^test_set_template>(
        test_set_function);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
