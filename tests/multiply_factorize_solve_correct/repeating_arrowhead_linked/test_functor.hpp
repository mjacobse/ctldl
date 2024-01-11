#include "tests/utility/generate_matrices.hpp"
#include "tests/utility/solution_generator.hpp"
#include "tests/utility/to_test_info.hpp"

#include <ctldl/factor_data/factorization_repeating_block_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_link.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_outer.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/matrix/multiply_repeating_block_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/matrix/matrix_link.hpp>
#include <ctldl/matrix/matrix_outer.hpp>
#include <ctldl/matrix/matrix_tridiagonal.hpp>
#include <ctldl/matrix/matrix_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/vector/block.hpp>
#include <ctldl/vector/flatten.hpp>
#include <ctldl/vector/vector_tridiagonal_arrowhead_linked.hpp>

#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace ctldl {

template <class TestMatrixTridiagDiag, class TestMatrixTridiagSubdiag,
          class PermutationTridiag, class TestMatrixLinkTridiag,
          class TestMatrixLinkDiag, class TestMatrixLinkOuter,
          class PermutationLink, class TestMatrixOuterSubdiag,
          class TestMatrixOuterDiag, class PermutationOuter, class Value,
          class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinked {
  void operator()(const SolutionGenerator& solution_generator,
                  const std::size_t num_repetitions,
                  std::mt19937& value_generator) const {
    constexpr auto dim_tridiag = std::size_t{TestMatrixTridiagDiag::Matrix::sparsity.num_rows};
    constexpr auto dim_link = std::size_t{TestMatrixLinkDiag::Matrix::sparsity.num_rows};
    constexpr auto dim_outer = std::size_t{TestMatrixOuterDiag::Matrix::sparsity.num_rows};

    constexpr Permutation<dim_tridiag> permutation_tridiag{PermutationTridiag::value};
    constexpr Permutation<dim_link> permutation_link{PermutationLink::value};
    constexpr Permutation<dim_outer> permutation_outer{PermutationOuter::value};

    constexpr auto sparsity = SparsityToFactorizeTridiagonalArrowheadLinked{
        SparsityToFactorizeTridiagonal{
            TestMatrixTridiagDiag::Matrix::sparsity,
            TestMatrixTridiagSubdiag::Matrix::sparsity,
            permutation_tridiag,
        },
        SparsityToFactorizeLink{
            TestMatrixLinkTridiag::Matrix::sparsity,
            TestMatrixLinkDiag::Matrix::sparsity,
            TestMatrixLinkOuter::Matrix::sparsity,
            permutation_link,
        },
        SparsityToFactorizeOuter{
            TestMatrixOuterSubdiag::Matrix::sparsity,
            TestMatrixOuterDiag::Matrix::sparsity,
            permutation_outer,
        },
    };

    FactorizationRepeatingBlockTridiagonalArrowheadLinked<sparsity, Value>
        factorization(num_repetitions);

    const MatrixTridiagonal matrices_tridiag{
        generateMatrices(TestMatrixTridiagDiag{}, value_generator, num_repetitions + 1),
        generateMatrices(TestMatrixTridiagSubdiag{}, value_generator, num_repetitions)};
    const MatrixLink matrices_link{
        generateMatrix(TestMatrixLinkTridiag{}, value_generator),
        generateMatrix(TestMatrixLinkDiag{}, value_generator),
        generateMatrix(TestMatrixLinkOuter{}, value_generator)};
    const MatrixOuter matrices_outer{
        generateMatrices(TestMatrixOuterSubdiag{}, value_generator, num_repetitions + 1),
        generateMatrix(TestMatrixOuterDiag{}, value_generator)};
    const MatrixTridiagonalArrowheadLinked matrix{
        matrices_tridiag, matrices_link, matrices_outer};
    factorization.factorize(matrix, FactorizeMethod{});

    const auto solution_tridiag = block<dim_tridiag, double>(
        solution_generator.generate((num_repetitions + 1) * dim_tridiag));
    const VectorTridiagonalArrowheadLinked solution{
        solution_tridiag, solution_generator.generate(dim_link),
        solution_generator.generate(dim_outer)};

    const std::vector<std::array<double, dim_tridiag>> rhs_tridiag(
        solution.tridiag.size());
    VectorTridiagonalArrowheadLinked rhs{rhs_tridiag,
                                         std::array<double, dim_link>{},
                                         std::array<double, dim_outer>{}};
    multiplyRepeatingBlockTridiagonalArrowheadLinked(matrix, solution, rhs);
    factorization.solveInPlace(rhs);

    BOOST_TEST_INFO_SCOPE("Matrix (tridiag, subdiag) " << TestMatrixTridiagDiag::description());
    BOOST_TEST_INFO_SCOPE("Matrix (tridiag, diag)    " << TestMatrixTridiagSubdiag::description());
    BOOST_TEST_INFO_SCOPE("Matrix (link, tridiag)    " << TestMatrixLinkTridiag::description());
    BOOST_TEST_INFO_SCOPE("Matrix (link, diag)       " << TestMatrixLinkDiag::description());
    BOOST_TEST_INFO_SCOPE("Matrix (link, outer)      " << TestMatrixLinkOuter::description());
    BOOST_TEST_INFO_SCOPE("Matrix (outer, subdiag):  " << TestMatrixOuterSubdiag::description());
    BOOST_TEST_INFO_SCOPE("Matrix (outer, diag):     " << TestMatrixOuterDiag::description());
    BOOST_TEST_INFO_SCOPE("Permutation (tridiag):    " << toTestInfo(permutation_tridiag));
    BOOST_TEST_INFO_SCOPE("Permutation (link):       " << toTestInfo(permutation_link));
    BOOST_TEST_INFO_SCOPE("Permutation (outer):      " << toTestInfo(permutation_outer));
    BOOST_TEST_INFO_SCOPE("Factor type:              " << toTestInfo(Value{}));
    BOOST_TEST_INFO_SCOPE("Num repetitions:          " << num_repetitions);
    BOOST_TEST_INFO_SCOPE("Solution:                 " << solution_generator.description);
    BOOST_TEST_INFO_SCOPE("Factorize method:         " << FactorizeMethod::description);

    const auto flatten = ctldl::flatten<double, dim_tridiag>;
    const auto tolerance =
        static_cast<double>(std::cbrt(std::numeric_limits<Value>::epsilon()));
    BOOST_TEST(flatten(rhs.tridiag) == flatten(solution.tridiag),
               boost::test_tools::tolerance(tolerance)
                   << boost::test_tools::per_element());
    BOOST_TEST(rhs.link == solution.link,
               boost::test_tools::tolerance(tolerance)
                   << boost::test_tools::per_element());
    BOOST_TEST(rhs.outer == solution.outer,
               boost::test_tools::tolerance(tolerance)
                   << boost::test_tools::per_element());
  }
};

}  // namespace ctldl
