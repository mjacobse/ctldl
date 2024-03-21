#include "tests/multiply_factorize_solve_correct/repeating_arrowhead_linked/test_functor_random.hpp"
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

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinked)

BOOST_AUTO_TEST_CASE(RandomExamples) {
  const auto seeds = TypeArgument<decltype(makeIntConstantSequence<8>())>{};

  constexpr int dim_start = 4;
  constexpr int dim_tridiag = 3;
  constexpr int dim_link = 2;
  constexpr int dim_outer = 5;
  const auto dimensions = makeTypeArgument<IntConstant<dim_start>>() ^
                          makeTypeArgument<IntConstant<dim_tridiag>>() ^
                          makeTypeArgument<IntConstant<dim_link>>() ^
                          makeTypeArgument<IntConstant<dim_outer>>();

  const auto factorize_value_types = makeTypeArgument<double>();
  const auto factorize_method =
      makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
  const auto solution_generators =
      makeValueArgument({getSolutionGeneratorNormallyDistributed(1000.0)});
  const auto repetition_counts = makeValueArgument<std::size_t>({0, 1, 2, 9});

  std::mt19937 value_generator{0};
  const auto test_set = seeds * dimensions * factorize_value_types *
                        factorize_method * solution_generators *
                        repetition_counts *
                        makeValueArgument({std::ref(value_generator)});
  foreach<TesterMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinkedRandom>(test_set);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
