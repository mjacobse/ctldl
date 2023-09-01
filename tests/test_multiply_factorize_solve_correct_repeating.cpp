#include "test_matrices/repeating/nos2.hpp"
#include "test_matrices/repeating/nos4.hpp"
#include "test_matrices/repeating/tridiagonal.hpp"
#include "utility/solution_generator.hpp"
#include "utility/test_set.hpp"
#include "utility/test_set_foreach.hpp"
#include "utility/to_test_info.hpp"

#include <ctldl/factor_data/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/matrix/multiply_repeating_block_tridiagonal.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/vector/block.hpp>
#include <ctldl/vector/flatten.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <tuple>
#include <vector>

namespace ctldl {
namespace {

template <class TestMatrix, class PermutationIn, class Value,
          class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrect {
  void operator()(const SolutionGenerator& solution_generator,
                  const std::size_t num_repetitions) const {
    constexpr auto& sparsity_A = TestMatrix::MatrixA::sparsity;
    constexpr auto& sparsity_B = TestMatrix::MatrixB::sparsity;

    constexpr auto block_dim = std::size_t{TestMatrix::block_dim};

    TestMatrix test_matrix(num_repetitions);
    FactorizationRepeatingBlockTridiagonal<sparsity_A, sparsity_B, Value,
                                           PermutationIn::permutation>
        factorization(num_repetitions);
    factorization.factorize(test_matrix.matrices_A, test_matrix.matrices_B,
                            FactorizeMethod{});

    const auto solution = block<block_dim, double>(
        solution_generator.generate((num_repetitions + 1) * block_dim));
    std::vector<std::array<double, block_dim>> rhs(solution.size());
    multiplyRepeatingBlockTridiagonal(test_matrix.matrices_A,
                                      test_matrix.matrices_B, solution, rhs);
    factorization.solveInPlace(rhs);

    const Permutation<block_dim> permutation(PermutationIn::permutation);
    BOOST_TEST_INFO("Test matrix:       " << test_matrix.description());
    BOOST_TEST_INFO("Permutation:       " << toTestInfo(permutation));
    BOOST_TEST_INFO("Factor type:       " << toTestInfo(Value{}));
    BOOST_TEST_INFO("Num repetitions:   " << num_repetitions);
    BOOST_TEST_INFO("Solution:          " << solution_generator.description);
    BOOST_TEST_INFO("Factorize method:  " << FactorizeMethod::description);

    const auto flatten = ctldl::flatten<double, block_dim>;
    const auto tolerance = TestMatrix::expected_error_amplifier *
                           std::numeric_limits<Value>::epsilon();
    BOOST_TEST(flatten(rhs) == flatten(solution),
               boost::test_tools::tolerance(tolerance)
                   << boost::test_tools::per_element());
  }
};

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectRepeating)

const auto solution_generators = makeValueArgument(
    {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
     getSolutionGeneratorNormallyDistributed(1000.0)});
const auto repetition_counts = makeValueArgument({0, 1, 2, 9});

const auto factorize_value_types = makeTypeArgument<double, float>();
const auto factorize_method =
    makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();

const auto common_test_set = factorize_value_types * factorize_method *
                             solution_generators * repetition_counts;

BOOST_AUTO_TEST_CASE(SmallExamplesAllPermutations) {
  const auto matrices_1x1 =
      makeTypeArgument<TestMatrixTridiagonal<1, double>>();
  const auto matrices_2x2 =
      makeTypeArgument<TestMatrixTridiagonal<2, double>>();
  const auto matrices_3x3 =
      makeTypeArgument<TestMatrixTridiagonal<3, double>,
                       TestMatrixTridiagonal<3, float>, TestMatrixNos2<double>,
                       TestMatrixNos2<float>>();
  const auto permutations_1x1 = TypeArgument<PermutationEnumeration<1>>{};
  const auto permutations_2x2 = TypeArgument<PermutationEnumeration<2>>{};
  const auto permutations_3x3 = TypeArgument<PermutationEnumeration<3>>{};
  const auto matrix_permutation_pairs = (matrices_1x1 * permutations_1x1) +
                                        (matrices_2x2 * permutations_2x2) +
                                        (matrices_3x3 * permutations_3x3);

  const auto test_set = matrix_permutation_pairs * common_test_set;
  foreach<TesterMultiplyFactorizeSolveCorrect>(test_set);
}

BOOST_AUTO_TEST_CASE(LargerExamplesGoodPermutation) {
  const auto matrix_permutation_pairs =
      makeTypeArgument<TestMatrixNos4<double>, TestMatrixNos4<float>>() *
      makeTypeArgument<TestPermutationNos4>();
  const auto test_set = matrix_permutation_pairs * common_test_set;
  foreach<TesterMultiplyFactorizeSolveCorrect>(test_set);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
