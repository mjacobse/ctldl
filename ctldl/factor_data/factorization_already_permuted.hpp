#pragma once

#include <ctldl/factor_data/factorization.hpp>
#include <ctldl/permutation/invert_permutation.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/sparsity/sparsity.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>

namespace ctldl {

// Basically Factorization, except that sparsity is already permuted according
// to permutation, so it should not be permuted again.
// We implement that by unpermuting the sparsity and letting the Factorization
// permute it back again. That way we factorize the correct sparsity and also
// remember the correct permutation within Factorization.
template <Sparsity sparsity, class Value, PermutationStatic permutation>
using FactorizationAlreadyPermuted =
    Factorization<getSparsityLowerTriangle<sparsity>(
                      invertPermutation(permutation)),
                  Value, permutation>;

}  // namespace ctldl
