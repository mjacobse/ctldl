#pragma once

#include <ctldl/permutation/permutation.hpp>

#include <string>

namespace ctldl {

std::string toTestInfo(PermutationView permutation);

inline const char* toTestInfo(double) { return "double"; }
inline const char* toTestInfo(float) { return "float"; }

}  // namespace ctldl
