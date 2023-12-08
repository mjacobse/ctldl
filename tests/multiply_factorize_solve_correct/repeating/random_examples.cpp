#include "tests/multiply_factorize_solve_correct/repeating/test_functor.hpp"
#include "tests/test_matrices/repeating/nos2.hpp"
#include "tests/test_matrices/repeating/tridiagonal.hpp"
#include "tests/utility/random/procedural_matrix_generator.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/test_set.hpp"
#include "tests/utility/test_set_foreach.hpp"

#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/permutation/permutation_identity.hpp>

#include <boost/test/unit_test.hpp>

#include <functional>
#include <random>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeating)

struct PermutationIdentityWrapped {
    static constexpr auto permutation = PermutationIdentity{};
};

BOOST_AUTO_TEST_CASE(RandomExamples) {
  const auto matrices_3x3 =
      TypeArgument<ProceduralMatrixGeneratorSymmetric<3>::Generate<0, 10>>{} *
      TypeArgument<ProceduralMatrixGenerator<3, 3>::Generate<1, 10>>{};
  const auto permutation_3x3 = makeTypeArgument<PermutationIdentityWrapped>();

  const auto factorize_value_types = makeTypeArgument<double>();
  const auto factorize_method =
      makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
  const auto solution_generators =
      makeValueArgument({getSolutionGeneratorNormallyDistributed(1000.0)});
  const auto repetition_counts = makeValueArgument<std::size_t>({0, 1, 2, 9});

  std::mt19937 value_generator{0};
  const auto test_set = matrices_3x3 * permutation_3x3 * factorize_value_types *
                        factorize_method * solution_generators *
                        repetition_counts *
                        makeValueArgument({std::ref(value_generator)});
  foreach<TesterMultiplyFactorizeSolveCorrectRepeating>(test_set);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
