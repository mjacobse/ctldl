#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace ctldl {

struct SolutionGenerator {
  std::function<std::vector<double>(std::size_t)> generate;
  const char* description;
};

SolutionGenerator getSolutionGeneratorAllOnes();
SolutionGenerator getSolutionGeneratorIota();
SolutionGenerator getSolutionGeneratorNormallyDistributed(
    double standard_deviation);

}  // namespace ctldl
