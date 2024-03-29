#include "tests/utility/random/procedural_test_permutation.hpp"
#include "tests/utility/test_static.hpp"

#include <ctldl/permutation/permutation.hpp>

#include <boost/test/unit_test.hpp>

#include <algorithm>

namespace ctldl {

BOOST_AUTO_TEST_CASE(ProceduralTestPermutationTest) {
  constexpr int seed = 0;
  constexpr int dim = 5;
  constexpr auto permutation = ProceduralTestPermutation<seed, dim>{}.value;
  constexpr Permutation<dim> permutation_identity;
  CTLDL_TEST_STATIC(std::is_permutation(
      permutation.cbegin(), permutation.cend(), permutation_identity.cbegin()));
}

}  // namespace ctldl
