#include "tests/utility/generate_matrices.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/to_test_info.hpp"

#include <ctldl/factor_data/factorization_repeating_block_tridiagonal.hpp>
#include <ctldl/matrix/multiply_repeating_block_tridiagonal.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/is_square.hpp>
#include <ctldl/vector/block.hpp>
#include <ctldl/vector/flatten.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace ctldl {

template <class TestMatrixA, class TestMatrixB,
          class PermutationIn, class Value, class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrectRepeating {
  void operator()(const SolutionGenerator& solution_generator,
                  const std::size_t num_repetitions,
                  std::mt19937& value_generator) const {
    constexpr auto& sparsity_A = TestMatrixA::Matrix::sparsity;
    constexpr auto& sparsity_B = TestMatrixB::Matrix::sparsity;

    static_assert(isSquare(sparsity_A));
    static_assert(isSquare(sparsity_B));
    static_assert(sparsity_B.num_rows == sparsity_A.num_rows);
    constexpr auto block_dim = std::size_t{sparsity_A.num_rows};
    constexpr Permutation<block_dim> permutation(PermutationIn::value);

    const auto matrices_A =
        generateMatrices(TestMatrixA{}, value_generator, num_repetitions + 1);
    const auto matrices_B =
        generateMatrices(TestMatrixB{}, value_generator, num_repetitions);
    FactorizationRepeatingBlockTridiagonal<sparsity_A, sparsity_B, Value,
                                           permutation>
        factorization(num_repetitions);
    factorization.factorize(matrices_A, matrices_B, FactorizeMethod{});

    const auto solution = block<block_dim, double>(
        solution_generator.generate((num_repetitions + 1) * block_dim));
    std::vector<std::array<double, block_dim>> rhs(solution.size());
    multiplyRepeatingBlockTridiagonal(matrices_A, matrices_B, solution, rhs);
    factorization.solveInPlace(rhs);

    BOOST_TEST_INFO("Test matrix A:     " << TestMatrixA::description());
    BOOST_TEST_INFO("Test matrix B:     " << TestMatrixB::description());
    BOOST_TEST_INFO("Permutation:       " << toTestInfo(permutation));
    BOOST_TEST_INFO("Factor type:       " << toTestInfo(Value{}));
    BOOST_TEST_INFO("Num repetitions:   " << num_repetitions);
    BOOST_TEST_INFO("Solution:          " << solution_generator.description);
    BOOST_TEST_INFO("Factorize method:  " << FactorizeMethod::description);

    const auto flatten = ctldl::flatten<double, block_dim>;
    const auto tolerance =
        static_cast<double>(std::cbrt(std::numeric_limits<Value>::epsilon()));
    BOOST_TEST(flatten(rhs) == flatten(solution),
               boost::test_tools::tolerance(tolerance)
                   << boost::test_tools::per_element());
  }
};

}  // namespace ctldl
