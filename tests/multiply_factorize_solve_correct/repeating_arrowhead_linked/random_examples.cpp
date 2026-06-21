#include "tests/multiply_factorize_solve_correct/repeating_arrowhead_linked/test_functor_random.hpp"
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

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinked)

BOOST_AUTO_TEST_CASE(RandomExamples) {
  static constexpr auto seeds = makeTemplateArgument(std::views::iota(0, 8));

  static constexpr int dim_start = 4;
  static constexpr int dim_tridiag = 3;
  static constexpr int dim_link = 2;
  static constexpr int dim_outer = 5;
  static constexpr auto dimensions =
      makeTemplateArgument(std::vector{dim_start}) ^
      makeTemplateArgument(std::vector{dim_tridiag}) ^
      makeTemplateArgument(std::vector{dim_link}) ^
      makeTemplateArgument(std::vector{dim_outer});

  static constexpr auto factorize_value_types = makeTemplateArgument<double>();
  static constexpr auto factorize_method =
      makeTemplateArgument<FactorizeMethodUpLooking,
                           FactorizeMethodEntryWise>();
  const auto solution_generators =
      makeFunctionArgument({getSolutionGeneratorNormallyDistributed(1000.0)});
  const auto repetition_counts =
      makeFunctionArgument<std::size_t>({0, 1, 2, 9});

  std::mt19937 value_generator{0};
  static constexpr auto test_set_template =
      seeds * dimensions * factorize_value_types * factorize_method;
  const auto test_set_function =
      solution_generators * repetition_counts *
      makeFunctionArgument({std::ref(value_generator)});
  foreach
    <^^TesterMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinkedRandom,
     ^^test_set_template>(test_set_function);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
