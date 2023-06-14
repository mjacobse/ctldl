#pragma once

#include <ctldl/permutation/permutation.hpp>

#include <sstream>
#include <string>

namespace ctldl {

template <std::size_t dim>
std::string toTestInfo(const Permutation<dim>& permutation) {
  static_assert(dim > 0);
  std::stringstream sstream;
  sstream << '[' << permutation[0];
  for (std::size_t i = 1; i < dim; ++i) {
    sstream << ' ' << permutation[i];
  }
  sstream << ']';
  return sstream.str();
}

const char* toTestInfo(double) { return "double"; }
const char* toTestInfo(float) { return "float"; }

}  // namespace ctldl