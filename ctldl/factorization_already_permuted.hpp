#pragma once

#include <ctldl/factorization.hpp>
#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>

namespace ctldl {

// Basically Factorization, except that sparsity is already permuted according
// to permutation, so it should not be permuted again.
// We implement that by unpermuting the sparsity and letting the Factorization
// permute it back again. That way we factorize the correct sparsity and also
// remember the correct permutation within Factorization.
template <auto sparsity, class Value, auto permutation>
using FactorizationAlreadyPermuted =
    Factorization<getSparsityLowerTriangle<sparsity>(
                      invertPermutation(permutation)),
                  Value, permutation>;

}  // namespace ctldl
