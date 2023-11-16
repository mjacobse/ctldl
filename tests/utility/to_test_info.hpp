#pragma once

#include <ctldl/permutation/permutation.hpp>

#include <sstream>
#include <string>

namespace ctldl {

template <std::size_t dim>
std::string toTestInfo(const Permutation<dim>& permutation) {
  if constexpr (dim == 0) {
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

inline const char* toTestInfo(double) { return "double"; }
inline const char* toTestInfo(float) { return "float"; }

}  // namespace ctldl
