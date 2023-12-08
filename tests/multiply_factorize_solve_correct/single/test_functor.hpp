#include "tests/utility/generate_matrices.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/to_test_info.hpp"

#include <ctldl/factor_data/factorization.hpp>
#include <ctldl/matrix/multiply.hpp>
#include <ctldl/permutation/permutation.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <limits>
#include <random>

namespace ctldl {

template <class TestMatrix, class PermutationIn, class Value,
          class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrectSingle {
  void operator()(const SolutionGenerator& solution_generator,
                  std::mt19937& value_generator) const {
    constexpr auto& sparsity = TestMatrix::Matrix::sparsity;
    constexpr auto dim = sparsity.num_rows;
    constexpr Permutation<dim> permutation(PermutationIn::permutation);

    auto test_matrix = generateMatrix(TestMatrix{}, value_generator);
    Factorization<sparsity, Value, permutation> factorization;
    factorization.factorize(test_matrix, FactorizeMethod{});

    const auto solution = solution_generator.generate(dim);
    std::array<double, dim> rhs = {0};
    multiply<MultiplyKind::Symmetric>(test_matrix, solution, rhs);
    factorization.solveInPlace(rhs);

    BOOST_TEST_INFO("Test matrix:       " << TestMatrix::description());
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

}  // namespace ctldl
