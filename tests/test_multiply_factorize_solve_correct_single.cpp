#include "test_matrices/repeating/nos2.hpp"
#include "test_matrices/repeating/tridiagonal.hpp"
#include "test_matrices/single/lfat5.hpp"
#include "utility/solution_generator.hpp"
#include "utility/test_set.hpp"
#include "utility/test_set_foreach.hpp"
#include "utility/to_test_info.hpp"

#include <ctldl/factor_data/factorization.hpp>
#include <ctldl/matrix/multiply.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/permutation/permutation_enumeration.hpp>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <limits>
#include <tuple>
#include <vector>

namespace ctldl {
namespace {

template <class TestMatrix, class PermutationIn, class Value,
          class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrect {
  void operator()(const SolutionGenerator& solution_generator) const {
    constexpr auto& sparsity = TestMatrix::sparsity;
    constexpr auto dim = sparsity.num_rows;
    constexpr Permutation<dim> permutation(PermutationIn::permutation);

    TestMatrix test_matrix;
    Factorization<sparsity, Value, permutation> factorization;
    factorization.factorize(test_matrix, FactorizeMethod{});

    const auto solution = solution_generator.generate(dim);
    std::array<double, dim> rhs = {0};
    multiply<MultiplyKind::Symmetric>(test_matrix, solution, rhs);
    factorization.solveInPlace(rhs);

    BOOST_TEST_INFO("Test matrix:       " << test_matrix.description());
    BOOST_TEST_INFO("Permutation:       " << toTestInfo(permutation));
    BOOST_TEST_INFO("Factor type:       " << toTestInfo(Value{}));
    BOOST_TEST_INFO("Solution:          " << solution_generator.description);
    BOOST_TEST_INFO("Factorize method:  " << FactorizeMethod::description);

    const auto tolerance = TestMatrix::expected_error_amplifier *
                           std::numeric_limits<Value>::epsilon();
    BOOST_TEST(rhs == solution, boost::test_tools::tolerance(tolerance)
                                    << boost::test_tools::per_element());
  }
};

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectSingle)

const auto solution_generators = makeValueArgument(
    {getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
     getSolutionGeneratorNormallyDistributed(1000.0)});

const auto factorize_value_types = makeTypeArgument<double, float>();
const auto factorize_method =
    makeTypeArgument<FactorizeMethodUpLooking, FactorizeMethodEntryWise>();
const auto common_test_set =
    factorize_value_types * factorize_method * solution_generators;

BOOST_AUTO_TEST_CASE(SmallExamplesAllPermutations) {
  const auto matrices_1x1 =
      makeTypeArgument<TestMatrixTridiagonal<1, double>::MatrixA>();
  const auto matrices_2x2 =
      makeTypeArgument<TestMatrixTridiagonal<2, double>::MatrixA>();
  const auto matrices_3x3 =
      makeTypeArgument<TestMatrixTridiagonal<3, double>::MatrixA,
                       TestMatrixTridiagonal<3, float>::MatrixA,
                       TestMatrixNos2<double>::MatrixA,
                       TestMatrixNos2<float>::MatrixA>();
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
      makeTypeArgument<TestMatrixLFAT5<double>, TestMatrixLFAT5<float>>() *
      makeTypeArgument<TestPermutationLFAT5>();

  const auto test_set = matrix_permutation_pairs * common_test_set;
  foreach<TesterMultiplyFactorizeSolveCorrect>(test_set);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
