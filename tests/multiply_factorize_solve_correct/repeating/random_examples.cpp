#include "tests/multiply_factorize_solve_correct/repeating/test_functor_random.hpp"
#include "tests/utility/int_constant.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/test_set.hpp"
#include "tests/utility/test_set_foreach.hpp"

#include <ctldl/factorize/factorize_method.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <functional>
#include <random>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeating)

BOOST_AUTO_TEST_CASE(RandomExamples) {
  static constexpr auto seeds =
      makeTypeArgumentFromTuple(makeIntConstantSequence<100>());
  static constexpr auto dimensions = makeTypeArgument<IntConstant<3>>();

  static constexpr auto factorize_value_types = makeTypeArgument<double>();
  static constexpr auto factorize_method =
      makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
  const auto solution_generators =
      makeValueArgument({getSolutionGeneratorNormallyDistributed(1000.0)});
  const auto repetition_counts = makeValueArgument<std::size_t>({0, 1, 2, 9});

  std::mt19937 value_generator{0};
  static constexpr auto test_set_types =
      seeds * dimensions * factorize_value_types * factorize_method;
  const auto test_set_values = solution_generators * repetition_counts *
                               makeValueArgument({std::ref(value_generator)});
  foreach
    <TesterMultiplyFactorizeSolveCorrectRepeatingRandom, ^^test_set_types>(
        test_set_values);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
