#include "solution_generator.hpp"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <random>
#include <vector>

namespace ctldl {

SolutionGenerator getSolutionGeneratorAllOnes() {
  const auto func = [](const std::size_t dim) {
    std::vector<double> ret(dim);
    std::fill(ret.begin(), ret.end(), 1.0);
    return ret;
  };
  return {func, "[1 1 ... 1 1]"};
}

SolutionGenerator getSolutionGeneratorIota() {
  const auto func = [](const std::size_t dim) {
    std::vector<double> ret(dim);
    std::iota(ret.begin(), ret.end(), 1.0);
    return ret;
  };
  return {func, "[1 2 ... d-1 d]"};
}

SolutionGenerator getSolutionGeneratorNormallyDistributed(
    const double standard_deviation) {
  const auto func = [standard_deviation](const std::size_t dim) {
    std::vector<double> ret(dim);
    const auto seed = dim;
    std::mt19937 rng(seed);
    std::normal_distribution<> distribution(0.0, standard_deviation);
    std::generate(ret.begin(), ret.end(),
                  [&distribution, &rng] { return distribution(rng); });
    return ret;
  };
  return {func, "normally distributed"};
}

}  // namespace ctldl
