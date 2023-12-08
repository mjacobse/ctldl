#include "tests/multiply_factorize_solve_correct/single/test_functor.hpp"
#include "tests/test_matrices/repeating/tridiagonal.hpp"
#include "tests/test_matrices/repeating/nos2.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/test_set.hpp"
#include "tests/utility/test_set_foreach.hpp"

#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/permutation/permutation_enumeration.hpp>

#include <boost/test/unit_test.hpp>

#include <functional>
#include <random>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectSingle)

BOOST_AUTO_TEST_CASE(SmallExamplesAllPermutations) {
  const auto matrices_1x1 =
      makeTypeArgument<TestMatrixTridiagonal<1, double>>();
  const auto matrices_2x2 =
      makeTypeArgument<TestMatrixTridiagonal<2, double>>();
  const auto matrices_3x3 =
      makeTypeArgument<TestMatrixTridiagonal<3, double>,
                       TestMatrixTridiagonal<3, float>,
                       TestMatrixNos2A<double>,
                       TestMatrixNos2A<float>>();
  const auto permutations_1x1 = TypeArgument<PermutationEnumeration<1>>{};
  const auto permutations_2x2 = TypeArgument<PermutationEnumeration<2>>{};
  const auto permutations_3x3 = TypeArgument<PermutationEnumeration<3>>{};
  const auto matrix_permutation_pairs = (matrices_1x1 * permutations_1x1) +
                                        (matrices_2x2 * permutations_2x2) +
                                        (matrices_3x3 * permutations_3x3);

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
