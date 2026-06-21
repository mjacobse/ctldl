#include "tests/multiply_factorize_solve_correct/repeating/test_functor.hpp"
#include "tests/test_matrices/repeating/indef.hpp"
#include "tests/test_matrices/repeating/nos2.hpp"
#include "tests/test_matrices/repeating/tridiagonal.hpp"
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

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeating)

BOOST_AUTO_TEST_CASE(SmallExamplesAllPermutations) {
  static constexpr auto matrices_1x1 =
      makeTemplateArgument<TestMatrixTridiagonal<1, double>>() *
      makeTemplateArgument<TestMatrixSingleEntryTopRight<1, double>>();
  static constexpr auto matrices_2x2 =
      makeTemplateArgument<TestMatrixTridiagonal<2, double>,
                           TestMatrixIndefA<double>>() ^
      makeTemplateArgument<TestMatrixSingleEntryTopRight<2, double>,
                           TestMatrixIndefB<double>>();
  static constexpr auto matrices_3x3 =
      makeTemplateArgument<TestMatrixTridiagonal<3, double>,
                           TestMatrixTridiagonal<3, float>,
                           TestMatrixNos2A<double>, TestMatrixNos2A<float>>() ^
      makeTemplateArgument<TestMatrixSingleEntryTopRight<3, double>,
                           TestMatrixSingleEntryTopRight<3, float>,
                           TestMatrixNos2B<double>, TestMatrixNos2B<float>>();
  static constexpr auto permutations_1x1 =
      makeTemplateArgument(enumeratePermutationsStatic(1));
  static constexpr auto permutations_2x2 =
      makeTemplateArgument(enumeratePermutationsStatic(2));
  static constexpr auto permutations_3x3 =
      makeTemplateArgument(enumeratePermutationsStatic(3));
  static constexpr auto matrix_permutation_pairs =
      (matrices_1x1 * permutations_1x1) + (matrices_2x2 * permutations_2x2) +
      (matrices_3x3 * permutations_3x3);

  static constexpr auto factorize_value_types =
      makeTemplateArgument<double, float>();
  static constexpr auto factorize_method =
      makeTemplateArgument<FactorizeMethodUpLooking,
                           FactorizeMethodEntryWise>();
  const auto solution_generators = makeFunctionArgument(
      {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
       getSolutionGeneratorNormallyDistributed(1000.0)});
  const auto repetition_counts =
      makeFunctionArgument<std::size_t>({0, 1, 2, 9});

  std::mt19937 value_generator{0};
  static constexpr auto test_set_template =
      matrix_permutation_pairs * factorize_value_types * factorize_method;
  const auto test_set_function =
      solution_generators * repetition_counts *
      makeFunctionArgument({std::ref(value_generator)});
  foreach
    <^^TesterMultiplyFactorizeSolveCorrectRepeating, ^^test_set_template>(
        test_set_function);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
