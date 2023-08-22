#include "test_matrices/repeating/nos2.hpp"
#include "test_matrices/repeating/tridiagonal.hpp"
#include "test_matrices/single/lfat5.hpp"
#include "utility/multi_invoke.hpp"
#include "utility/solution_generator.hpp"
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
struct TesterMultiplyFactorizeSolveCorrectImpl {
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

template <class TestCaseFactorize, class Value, class FactorizeMethod>
using TesterMultiplyFactorizeSolveCorrect =
    TesterMultiplyFactorizeSolveCorrectImpl<
        typename TestCaseFactorize::Matrix,
        typename TestCaseFactorize::Permutation, Value, FactorizeMethod>;

template <class TestMatrix, class TestPermutation>
struct TestCaseFactorize {
  using Matrix = TestMatrix;
  using Permutation = TestPermutation;
};

template <class TestMatrix, class... TestPermutations>
auto getTestCasesFactorize(std::tuple<TestPermutations...>) {
  return std::tuple<TestCaseFactorize<TestMatrix, TestPermutations>...>{};
}

template <class TestMatrix>
auto getTestCasesFactorizeAllPermutations() {
  using Permutations = PermutationEnumeration<TestMatrix::sparsity.num_rows>;
  return getTestCasesFactorize<TestMatrix>(Permutations{});
}

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrectSingle)

BOOST_AUTO_TEST_CASE(SmallExamplesAllPermutations) {
  using TestCasesFactorize = decltype(std::tuple_cat(
      getTestCasesFactorizeAllPermutations<TestMatrixNos2<double>::MatrixA>(),
      getTestCasesFactorizeAllPermutations<TestMatrixNos2<float>::MatrixA>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<1, double>::MatrixA>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<2, double>::MatrixA>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<3, double>::MatrixA>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<3, float>::MatrixA>()));
  using FactorizeValueTypes = std::tuple<double, float>;
  using FactorizeMethods =
      std::tuple<FactorizeMethodUpLooking, FactorizeMethodEntryWise>;
  const std::vector<SolutionGenerator> solution_generators = {
      getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
      getSolutionGeneratorNormallyDistributed(1000.0)};

  using TestTypes =
      std::tuple<TestCasesFactorize, FactorizeValueTypes, FactorizeMethods>;
  const auto test_values = std::make_tuple(solution_generators);
  multiInvoke<TesterMultiplyFactorizeSolveCorrect, TestTypes>(test_values);
}

BOOST_AUTO_TEST_CASE(LargerExamplesGoodPermutation) {
  using TestCasesFactorize = std::tuple<
      TestCaseFactorize<TestMatrixLFAT5<double>, TestPermutationLFAT5>,
      TestCaseFactorize<TestMatrixLFAT5<float>, TestPermutationLFAT5>>;
  using FactorizeValueTypes = std::tuple<double, float>;
  using FactorizeMethods =
      std::tuple<FactorizeMethodUpLooking, FactorizeMethodEntryWise>;
  const std::vector<SolutionGenerator> solution_generators = {
      getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
      getSolutionGeneratorNormallyDistributed(1000.0)};

  using TestTypes =
      std::tuple<TestCasesFactorize, FactorizeValueTypes, FactorizeMethods>;
  const auto test_values = std::make_tuple(solution_generators);
  multiInvoke<TesterMultiplyFactorizeSolveCorrect, TestTypes>(test_values);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
