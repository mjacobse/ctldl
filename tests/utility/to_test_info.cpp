#include "to_test_info.hpp"

#include <sstream>

namespace ctldl {

std::string toTestInfo(const PermutationView permutation) {
  const auto dim = std::size_t{permutation.size()};
  if (dim == 0) {
    return "[]";
  }

  std::stringstream sstream;
  sstream << '[' << permutation[0];
  for (std::size_t i = 1; i < dim; ++i) {
    sstream << ' ' << permutation[i];
  }
  sstream << ']';
  return sstream.str();
}

}  // namespace ctldl
