#pragma once

#include "tests/multiply_factorize_solve_correct/repeating_arrowhead_linked/test_functor.hpp"
#include "tests/utility/random/procedural_test_matrix.hpp"
#include "tests/utility/random/procedural_test_permutation.hpp"
#include "tests/utility/solution_generator.hpp"

#include <cstddef>
#include <random>

namespace ctldl {

template <class Seed, class DimTridiag, class DimLink, class DimOuter,
          class Value, class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinkedRandom {
  void operator()(const SolutionGenerator& solution_generator,
                  const std::size_t num_repetitions,
                  std::mt19937& value_generator) const {
    constexpr auto dim_tridiag = DimTridiag::value;
    constexpr auto dim_link = DimLink::value;
    constexpr auto dim_outer = DimOuter::value;

    constexpr auto seed0 = Seed::value;
    using MatrixTridiagDiag = ProceduralTestMatrixSymmetricT<seed0, dim_tridiag>;
    constexpr auto seed1 = MatrixTridiagDiag::next_generator.state;
    using MatrixTridiagSubdiag = ProceduralTestMatrixT<seed1, dim_tridiag, dim_tridiag>;
    constexpr auto seed2 = MatrixTridiagSubdiag::next_generator.state;
    using MatrixLinkTridiag = ProceduralTestMatrixT<seed2, dim_link, dim_tridiag>;
    constexpr auto seed3 = MatrixLinkTridiag::next_generator.state;
    using MatrixLinkDiag = ProceduralTestMatrixSymmetricT<seed3, dim_link>;
    constexpr auto seed4 = MatrixLinkDiag::next_generator.state;
    using MatrixLinkOuter = ProceduralTestMatrixT<seed4, dim_outer, dim_link>;
    constexpr auto seed5 = MatrixLinkOuter::next_generator.state;
    using MatrixOuterSubdiag = ProceduralTestMatrixT<seed5, dim_outer, dim_tridiag>;
    constexpr auto seed6 = MatrixOuterSubdiag::next_generator.state;
    using MatrixOuterDiag = ProceduralTestMatrixSymmetricT<seed6, dim_outer>;

    constexpr auto seed7 = MatrixOuterDiag::next_generator.state;
    using PermutationTridiag = ProceduralTestPermutation<seed7, dim_tridiag>;
    constexpr auto seed8 = PermutationTridiag::next_generator.state;
    using PermutationLink = ProceduralTestPermutation<seed8, dim_link>;
    constexpr auto seed9 = PermutationLink::next_generator.state;
    using PermutationOuter = ProceduralTestPermutation<seed9, dim_outer>;

    TesterMultiplyFactorizeSolveCorrectRepeatingBlockTridiagonalArrowheadLinked<
        MatrixTridiagDiag, MatrixTridiagSubdiag, PermutationTridiag,
        MatrixLinkTridiag, MatrixLinkDiag, MatrixLinkOuter, PermutationLink,
        MatrixOuterSubdiag, MatrixOuterDiag, PermutationOuter, Value,
        FactorizeMethod>{}(solution_generator, num_repetitions,
                           value_generator);
  }
};

}  // namespace ctldl
