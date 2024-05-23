#include <ctldl/factor_data/factorization_repeating_block_tridiagonal_arrowhead_linked.hpp>

#include "tests/utility/test_static.hpp"

#include <boost/test/unit_test.hpp>

#include <array>
#include <cstddef>

namespace ctldl {
namespace {

BOOST_AUTO_TEST_SUITE(FactorSparsity)

BOOST_AUTO_TEST_CASE(SplinePrepermuted) {
  constexpr auto dim_start = std::size_t{4};
  constexpr auto dim_tridiag = std::size_t{7};
  constexpr auto dim_link = std::size_t{0};
  constexpr auto dim_outer = std::size_t{0};

  constexpr auto sparsity = SparsityToFactorizeTridiagonalArrowheadLinked{
      SparsityToFactorizeStart{
          makeEmptySparsity<dim_start, dim_start>(),
          makeSparsity<dim_tridiag, dim_start>(
              {{3, 0}, {3, 1}, {4, 1}, {4, 3}, {5, 2}, {5, 3}}),
          makeEmptySparsity<dim_outer, dim_start>()},
      SparsityToFactorizeTridiagonal{
          makeSparsity<dim_tridiag, dim_tridiag>(
              {{3, 0}, {3, 2}, {4, 2}, {5, 1}, {6, 4}, {6, 5}}),
          makeSparsity<dim_tridiag, dim_tridiag>(
              {{3, 0}, {3, 2}, {4, 2}, {4, 6}, {5, 1}, {5, 6}})},
      SparsityToFactorizeLink{makeEmptySparsity<dim_link, dim_tridiag>(),
                              makeEmptySparsity<dim_link, dim_link>(),
                              makeEmptySparsity<dim_outer, dim_link>()},
      SparsityToFactorizeOuter{makeEmptySparsity<dim_outer, dim_tridiag>(),
                               makeEmptySparsity<dim_outer, dim_outer>()}};

  using Factorization =
      FactorizationRepeatingBlockTridiagonalArrowheadLinked<sparsity, double>;

  constexpr auto diag_correct = makeSparsity<dim_tridiag, dim_tridiag>(
      {{3, 0}, {3, 2}, {4, 2}, {4, 3}, {5, 1}, {5, 3}, {5, 4}, {6, 4}, {6, 5}});
  CTLDL_TEST_STATIC(isSparsityEqual(Factorization::FactorTridiagDiag::sparsity,
                                    diag_correct));
  constexpr auto subdiag_correct = makeSparsity<dim_tridiag, dim_tridiag>(
    {{3, 0},         {3, 2}, {3, 3}, {3, 4}, {3, 5}, {3, 6},
                     {4, 2}, {4, 3}, {4, 4}, {4, 5}, {4, 6},
             {5, 1},                         {5, 5}, {5, 6}});
  CTLDL_TEST_STATIC(isSparsityEqual(
      Factorization::FactorTridiagSubdiag::sparsity, subdiag_correct));
}

BOOST_AUTO_TEST_CASE(Spline) {
  constexpr auto dim_start = std::size_t{4};
  constexpr auto dim_tridiag = std::size_t{7};
  constexpr auto dim_link = std::size_t{0};
  constexpr auto dim_outer = std::size_t{0};

  constexpr Permutation permutation_start{
      std::array<std::size_t, dim_start>{0, 1, 2, 3}};
  constexpr Permutation permutation_tridiag{
      std::array<std::size_t, dim_tridiag>{3, 5, 4, 0, 1, 2, 6}};

  constexpr auto sparsity = SparsityToFactorizeTridiagonalArrowheadLinked{
      SparsityToFactorizeStart{
          makeEmptySparsity<dim_start, dim_start>(),
          makeSparsity<dim_tridiag, dim_start>(
              {{0, 0}, {0, 1}, {1, 1}, {1, 3}, {2, 2}, {2, 3}}),
          makeEmptySparsity<dim_outer, dim_start>(),
          permutation_start,
      },
      SparsityToFactorizeTridiagonal{
          makeSparsity<dim_tridiag, dim_tridiag>(
              {{3, 0}, {4, 0}, {4, 1}, {5, 2}, {6, 1}, {6, 2}}),
          makeSparsity<dim_tridiag, dim_tridiag>(
              {{0, 3}, {0, 4}, {1, 4}, {1, 6}, {2, 5}, {2, 6}}),
          permutation_tridiag,
      },
      SparsityToFactorizeLink{makeEmptySparsity<dim_link, dim_tridiag>(),
                              makeEmptySparsity<dim_link, dim_link>(),
                              makeEmptySparsity<dim_outer, dim_link>()},
      SparsityToFactorizeOuter{makeEmptySparsity<dim_outer, dim_tridiag>(),
                               makeEmptySparsity<dim_outer, dim_outer>()}};

  using Factorization =
      FactorizationRepeatingBlockTridiagonalArrowheadLinked<sparsity, double>;

  constexpr auto diag_correct = makeSparsity<dim_tridiag, dim_tridiag>(
      {{3, 0}, {3, 2}, {4, 2}, {4, 3}, {5, 1}, {5, 3}, {5, 4}, {6, 4}, {6, 5}});
  CTLDL_TEST_STATIC(isSparsityEqual(Factorization::FactorTridiagDiag::sparsity,
                                    diag_correct));
  constexpr auto subdiag_correct = makeSparsity<dim_tridiag, dim_tridiag>(
    {{3, 0},         {3, 2}, {3, 3}, {3, 4}, {3, 5}, {3, 6},
                     {4, 2}, {4, 3}, {4, 4}, {4, 5}, {4, 6},
             {5, 1},                         {5, 5}, {5, 6}});
  CTLDL_TEST_STATIC(isSparsityEqual(
      Factorization::FactorTridiagSubdiag::sparsity, subdiag_correct));
}

BOOST_AUTO_TEST_SUITE_END()

}  // anonymous namespace
}  // namespace ctldl
