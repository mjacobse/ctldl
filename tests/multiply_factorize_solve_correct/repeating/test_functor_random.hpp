#pragma once

#include "tests/multiply_factorize_solve_correct/repeating/test_functor.hpp"
#include "tests/utility/random/procedural_test_matrix.hpp"
#include "tests/utility/random/procedural_test_permutation.hpp"
#include "tests/utility/solution_generator.hpp"

#include <cstddef>
#include <random>

namespace ctldl {

template <class Seed, class Dim, class Value, class FactorizeMethod>
struct TesterMultiplyFactorizeSolveCorrectRepeatingRandom {
  void operator()(const SolutionGenerator& solution_generator,
                  const std::size_t num_repetitions,
                  std::mt19937& value_generator) const {
    constexpr auto dim = Dim::value;

    constexpr auto seed0 = Seed::value;
    using MatrixDiag = ProceduralTestMatrixSymmetricT<seed0, dim>;
    constexpr auto seed1 = MatrixDiag::next_generator.state;
    using MatrixSubdiag = ProceduralTestMatrixT<seed1, dim, dim>;

    constexpr auto seed2 = MatrixSubdiag::next_generator.state;
    using PermutationTridiag = ProceduralTestPermutation<seed2, dim>;

    TesterMultiplyFactorizeSolveCorrectRepeating<MatrixDiag, MatrixSubdiag,
                                                 PermutationTridiag, Value,
                                                 FactorizeMethod>{}(
        solution_generator, num_repetitions, value_generator);
  }
};

}  // namespace ctldl
