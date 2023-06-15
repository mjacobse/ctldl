#include "test_matrices/nos2.hpp"
#include "test_matrices/nos4.hpp"
#include "test_matrices/tridiagonal.hpp"
#include "utility/multi_invoke.hpp"
#include "utility/solution_generator.hpp"
#include "utility/to_test_info.hpp"

#include <ctldl/factorization_repeating_block_tridiagonal.hpp>
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
struct TesterMultiplyFactorizeSolveCorrectImpl {
  void operator()(const SolutionGenerator& solution_generator,
                  const std::size_t num_repetitions) const {
    using SparsityA = typename TestMatrix::MatrixA::Sparsity;
    using SparsityB = typename TestMatrix::MatrixB::Sparsity;

    constexpr auto block_dim = std::size_t{TestMatrix::dim};

    TestMatrix test_matrix(num_repetitions);
    FactorizationRepeatingBlockTridiagonal<SparsityA, SparsityB, Value,
                                           PermutationIn>
        factorization(num_repetitions);
    factorization.factor(test_matrix.matrices_A, test_matrix.matrices_B,
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
  using Permutations = PermutationEnumeration<TestMatrix::dim>;
  return getTestCasesFactorize<TestMatrix>(Permutations{});
}

BOOST_AUTO_TEST_SUITE(TestMultiplyFactorizeSolveCorrect)

BOOST_AUTO_TEST_CASE(SmallExamplesAllPermutations) {
  using TestCasesFactorize = decltype(std::tuple_cat(
      getTestCasesFactorizeAllPermutations<TestMatrixNos2<double>>(),
      getTestCasesFactorizeAllPermutations<TestMatrixNos2<float>>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<1, double>>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<2, double>>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<3, double>>(),
      getTestCasesFactorizeAllPermutations<TestMatrixTridiagonal<3, float>>()));
  using FactorizeValueTypes = std::tuple<double, float>;
  using FactorizeMethods =
      std::tuple<FactorizeMethodUpLooking, FactorizeMethodEntryWise>;
  const std::vector<SolutionGenerator> solution_generators = {
      getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
      getSolutionGeneratorNormallyDistributed(1000.0)};
  const std::vector<std::size_t> repetition_counts = {0, 1, 2, 9};

  using TestTypes =
      std::tuple<TestCasesFactorize, FactorizeValueTypes, FactorizeMethods>;
  const auto test_values =
      std::make_tuple(solution_generators, repetition_counts);
  multiInvoke<TesterMultiplyFactorizeSolveCorrect, TestTypes>(test_values);
}

BOOST_AUTO_TEST_CASE(LargerExamplesGoodPermutation) {
  using TestCasesFactorize =
      std::tuple<TestCaseFactorize<TestMatrixNos4<double>, TestPermutationNos4>,
                 TestCaseFactorize<TestMatrixNos4<float>, TestPermutationNos4>>;
  using FactorizeValueTypes = std::tuple<double, float>;
  using FactorizeMethods =
      std::tuple<FactorizeMethodUpLooking, FactorizeMethodEntryWise>;
  const std::vector<SolutionGenerator> solution_generators = {
      getSolutionGeneratorAllOnes(), getSolutionGeneratorIota(),
      getSolutionGeneratorNormallyDistributed(1000.0)};
  const std::vector<std::size_t> repetition_counts = {0, 1, 2, 9};

  using TestTypes =
      std::tuple<TestCasesFactorize, FactorizeValueTypes, FactorizeMethods>;
  const auto test_values =
      std::make_tuple(solution_generators, repetition_counts);
  multiInvoke<TesterMultiplyFactorizeSolveCorrect, TestTypes>(test_values);
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
