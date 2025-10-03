#pragma once

#include "tests/utility/random/compile_time_generate.hpp"
#include "tests/utility/random/linear_congruence_generator.hpp"
#include "tests/utility/random/shuffle.hpp"

#include <ctldl/permutation/permutation.hpp>

#include <array>
#include <cstddef>

namespace ctldl {

template <std::size_t dim>
struct PermutationDistribution {
  using result_type = PermutationStatic<dim>;

  template <class Generator>
  constexpr auto operator()(Generator& generator) const {
    std::array<std::size_t, dim> indices = PermutationStatic<dim>().indices();
    shuffle(indices.begin(), indices.end(), generator);
    return PermutationStatic<dim>(indices);
  }
};

template <std::size_t seed, std::size_t dim>
struct ProceduralTestPermutation {
 private:
  static constexpr auto permutation_generated = compileTimeGenerate(
      PermutationDistribution<dim>{}, LinearCongruenceGenerator{seed});

 public:
  static constexpr auto value = permutation_generated.result;
  static constexpr auto next_generator = permutation_generated.generator;
};

}  // namespace ctldl
