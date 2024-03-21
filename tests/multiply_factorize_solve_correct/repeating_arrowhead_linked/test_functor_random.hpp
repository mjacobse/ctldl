#pragma once

#include "tests/multiply_factorize_solve_correct/repeating_arrowhead_linked/test_functor.hpp"
#include "tests/utility/random/procedural_test_matrix.hpp"
#include "tests/utility/random/procedural_test_permutation.hpp"
#include "tests/utility/solution_generator.hpp"

#include <cstddef>
#include <random>

namespace ctldl {

template <class Seed, class DimStart, class DimTridiag, class DimLink,
          class DimOuter, class Value, class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinkedRandom {
  void operator()(const SolutionGenerator& solution_generator,
                  const std::size_t num_repetitions,
                  std::mt19937& value_generator) const {
    constexpr auto dim_start = DimStart::value;
    constexpr auto dim_tridiag = DimTridiag::value;
    constexpr auto dim_link = DimLink::value;
    constexpr auto dim_outer = DimOuter::value;

    constexpr auto seed0 = Seed::value;
    using MatrixStartDiag = ProceduralTestMatrixSymmetricT<seed0, dim_start>;
    constexpr auto seed1 = MatrixStartDiag::next_generator.state;
    using MatrixStartTridiag = ProceduralTestMatrixT<seed1, dim_tridiag, dim_start>;
    constexpr auto seed2 = MatrixStartTridiag::next_generator.state;
    using MatrixStartOuter = ProceduralTestMatrixT<seed2, dim_outer, dim_start>;
    constexpr auto seed3 = MatrixStartOuter::next_generator.state;
    using MatrixTridiagDiag = ProceduralTestMatrixSymmetricT<seed3, dim_tridiag>;
    constexpr auto seed4 = MatrixTridiagDiag::next_generator.state;
    using MatrixTridiagSubdiag = ProceduralTestMatrixT<seed4, dim_tridiag, dim_tridiag>;
    constexpr auto seed5 = MatrixTridiagSubdiag::next_generator.state;
    using MatrixLinkTridiag = ProceduralTestMatrixT<seed5, dim_link, dim_tridiag>;
    constexpr auto seed6 = MatrixLinkTridiag::next_generator.state;
    using MatrixLinkDiag = ProceduralTestMatrixSymmetricT<seed6, dim_link>;
    constexpr auto seed7 = MatrixLinkDiag::next_generator.state;
    using MatrixLinkOuter = ProceduralTestMatrixT<seed7, dim_outer, dim_link>;
    constexpr auto seed8 = MatrixLinkOuter::next_generator.state;
    using MatrixOuterSubdiag = ProceduralTestMatrixT<seed8, dim_outer, dim_tridiag>;
    constexpr auto seed9 = MatrixOuterSubdiag::next_generator.state;
    using MatrixOuterDiag = ProceduralTestMatrixSymmetricT<seed9, dim_outer>;

    constexpr auto seed10 = MatrixOuterDiag::next_generator.state;
    using PermutationStart = ProceduralTestPermutation<seed10, dim_start>;
    constexpr auto seed11 = MatrixOuterDiag::next_generator.state;
    using PermutationTridiag = ProceduralTestPermutation<seed11, dim_tridiag>;
    constexpr auto seed12 = PermutationTridiag::next_generator.state;
    using PermutationLink = ProceduralTestPermutation<seed12, dim_link>;
    constexpr auto seed13 = PermutationLink::next_generator.state;
    using PermutationOuter = ProceduralTestPermutation<seed13, dim_outer>;

    TesterMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinked<
        MatrixStartDiag, MatrixStartTridiag, MatrixStartOuter, PermutationStart,
        MatrixTridiagDiag, MatrixTridiagSubdiag, PermutationTridiag,
        MatrixLinkTridiag, MatrixLinkDiag, MatrixLinkOuter, PermutationLink,
        MatrixOuterSubdiag, MatrixOuterDiag, PermutationOuter, Value,
        FactorizeMethod>{}(solution_generator, num_repetitions,
                           value_generator);
  }
};

}  // namespace ctldl
