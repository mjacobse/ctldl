#pragma once

#include <ctldl/factor_data/factorization_already_permuted.hpp>
#include <ctldl/factor_data/factorization_subdiagonal_block.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factorize/factorize_entry_wise.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/factorize/factorize_up_looking.hpp>
#include <ctldl/factorize/regularization.hpp>
#include <ctldl/matrix/matrix_link.hpp>
#include <ctldl/matrix/matrix_outer.hpp>
#include <ctldl/matrix/matrix_start.hpp>
#include <ctldl/matrix/matrix_tridiagonal.hpp>
#include <ctldl/matrix/matrix_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/solve/solve_backward_substitution.hpp>
#include <ctldl/solve/solve_forward_substitution.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/filled_in_sparsity_blocked.hpp>
#include <ctldl/symbolic/filled_in_sparsity_repeating.hpp>
#include <ctldl/utility/contracts.hpp>

#include <cstddef>
#include <limits>
#include <source_location>
#include <span>
#include <vector>


namespace ctldl {

/**
 * Factorization of a symmetric block-tridiagonal arrowhead matrix with
 * repeating sparsity and connecting link blocks between the tridiagonal and
 * outer part. See \ref SparsityToFactorizeTridiagonalArrowheadLinked
 */
template <SparsityToFactorizeTridiagonalArrowheadLinked sparsity,
          class Value_>
class FactorizationRepeatingBlockTridiagonalArrowheadLinked {
 private:
  static constexpr auto dim_start = std::size_t{sparsity.dim_start};
  static constexpr auto dim_tridiag = std::size_t{sparsity.dim_tridiag};
  static constexpr auto dim_link = std::size_t{sparsity.dim_link};
  static constexpr auto dim_outer = std::size_t{sparsity.dim_outer};
  static constexpr auto permutation_start = sparsity.start.permutation;
  static constexpr auto permutation_tridiag = sparsity.tridiag.permutation;
  static constexpr auto permutation_link = sparsity.link.permutation;
  static constexpr auto permutation_outer = sparsity.outer.permutation;

  static constexpr auto helper_start = getFilledInSparsityBlocked3x3<
      sparsity.start.diag, sparsity.start.next, sparsity.tridiag.diag,
      sparsity.start.outer, sparsity.outer.subdiag, sparsity.outer.diag,
      permutation_start, permutation_tridiag, permutation_outer>();
  static constexpr auto sparsity_factor_start_diag = helper_start.block11;
  static constexpr auto sparsity_factor_start_tridiag = helper_start.block21;
  static constexpr auto sparsity_factor_start_outer = helper_start.block31;
  static constexpr auto sparsity_tridiag_diag = helper_start.block22;
  static constexpr auto sparsity_tridiag_subdiag = getSparsityStaticPermuted(
      sparsity.tridiag.subdiag, permutation_tridiag, permutation_tridiag);
  static constexpr auto sparsity_outer_subdiag = helper_start.block32;
  static constexpr auto sparsity_outer_diag = helper_start.block33;

  static constexpr auto sparsity_factor = getFilledInSparsityRepeatingArrowhead<
      sparsity_tridiag_diag, sparsity_tridiag_subdiag, sparsity_outer_subdiag,
      PermutationStatic<dim_tridiag>{}, PermutationStatic<dim_outer>{}>();

  static constexpr auto helper = [] {
    constexpr auto sparsity_link_tridiag_permuted_cols =
        getSparsityStaticPermuted(sparsity.link.prev,
                                  PermutationStatic<dim_link>{},
                                  permutation_tridiag);
    constexpr auto sparsity_link_outer_permuted_rows =
        getSparsityStaticPermuted(sparsity.link.next, permutation_outer,
                                  PermutationStatic<dim_link>{});
    return getFilledInSparsityBlocked3x3<
        sparsity_factor.diag, sparsity_link_tridiag_permuted_cols,
        sparsity.link.diag, sparsity_factor.outer,
        sparsity_link_outer_permuted_rows, sparsity_outer_diag,
        PermutationStatic<dim_tridiag>{}, permutation_link,
        PermutationStatic<dim_outer>{}>();
  }();
  static constexpr auto sparsity_factor_link_tridiag = helper.block21;
  static constexpr auto sparsity_factor_link_diag = helper.block22;
  static constexpr auto sparsity_factor_link_outer = helper.block32;
  static constexpr auto sparsity_factor_outer_diag = helper.block33;

 public:
  using Value = Value_;
  using FactorStartDiag =
      FactorizationAlreadyPermuted<sparsity_factor_start_diag, Value,
                                   permutation_start>;
  using FactorStartTridiag =
      FactorizationSubdiagonalBlock<sparsity_factor_start_tridiag, Value,
                                    permutation_tridiag, permutation_start>;
  using FactorStartOuter =
      FactorizationSubdiagonalBlock<sparsity_factor_start_outer, Value,
                                    permutation_outer, permutation_start>;
  using FactorTridiagDiag =
      FactorizationAlreadyPermuted<sparsity_factor.diag, Value,
                                   permutation_tridiag>;
  using FactorTridiagSubdiag =
      FactorizationSubdiagonalBlock<sparsity_factor.subdiag, Value,
                                    permutation_tridiag, permutation_tridiag>;
  using FactorLinkTridiag =
      FactorizationSubdiagonalBlock<sparsity_factor_link_tridiag, Value,
                                    permutation_link, permutation_tridiag>;
  using FactorLinkDiag = FactorizationAlreadyPermuted<sparsity_factor_link_diag,
                                                      Value, permutation_link>;
  using FactorLinkOuter =
      FactorizationSubdiagonalBlock<sparsity_factor_link_outer, Value,
                                    permutation_outer, permutation_link>;
  using FactorOuterSubdiag =
      FactorizationSubdiagonalBlock<sparsity_factor.outer, Value,
                                    permutation_outer, permutation_tridiag>;
  using FactorOuterDiag =
      FactorizationAlreadyPermuted<sparsity_factor_outer_diag, Value,
                                   permutation_outer>;

  using FactorStart =
      MatrixStart<FactorStartDiag, FactorStartTridiag, FactorStartOuter>;
  using FactorTridiagonal =
      MatrixTridiagonal<std::vector<FactorTridiagDiag>,
                        std::vector<FactorTridiagSubdiag>>;
  using FactorLink =
      MatrixLink<FactorLinkTridiag, FactorLinkDiag, FactorLinkOuter>;
  using FactorOuter =
      MatrixOuter<std::vector<FactorOuterSubdiag>, FactorOuterDiag>;

  explicit FactorizationRepeatingBlockTridiagonalArrowheadLinked(
      const std::size_t num_repetitions)
      : m_data(MatrixTridiagonalArrowheadLinked{
            FactorStart{},
            MatrixTridiagonal{
                std::vector<FactorTridiagDiag>(num_repetitions + 1),
                std::vector<FactorTridiagSubdiag>(num_repetitions)},
            FactorLink{},
            MatrixOuter{std::vector<FactorOuterSubdiag>(num_repetitions + 1),
                        FactorOuterDiag{}}}) {
    checkNumRepetitionsClassInvariant();
  }

  void checkNumRepetitionsClassInvariant(
      std::source_location source_location =
          std::source_location::current()) const {
    contract_assert(numRepetitions() < std::numeric_limits<std::size_t>::max(),
                    source_location);
    contract_assert(m_data.tridiag.diag.size() == numRepetitions() + 1,
                    source_location);
    contract_assert(m_data.tridiag.subdiag.size() == numRepetitions(),
                    source_location);
    contract_assert(m_data.outer.subdiag.size() == numRepetitions() + 1,
                    source_location);
  }

  auto numRepetitions() const noexcept { return m_data.tridiag.subdiag.size(); }
  const auto& data() const noexcept { return m_data; };

  template <class FactorizeMethodTag = FactorizeMethodUpLooking,
            class MatrixBlocksInputValues>
  void factorize(const MatrixBlocksInputValues& input,
                 const Regularization auto& regularization,
                 const FactorizeMethodTag method_tag = {}) {
    const auto num_repetitions = std::size_t{numRepetitions()};
    // start
    m_data.start.diag.factorize(input.start.diag, regularization, method_tag);
    // tridiagonal
    ::ctldl::factorize(m_data.start.diag, input.start.next,
                       input.tridiag.diag[0], m_data.start.next,
                       m_data.tridiag.diag[0], regularization, method_tag);
    for (std::size_t i = 0; i < num_repetitions; ++i) {
      ::ctldl::factorize(m_data.tridiag.diag[i], input.tridiag.subdiag[i],
                         input.tridiag.diag[i + 1], m_data.tridiag.subdiag[i],
                         m_data.tridiag.diag[i + 1], regularization,
                         method_tag);
    }
    ::ctldl::factorize(m_data.tridiag.diag[num_repetitions], input.link.prev,
                       input.link.diag, m_data.link.prev, m_data.link.diag,
                       regularization, method_tag);
    // outer
    const FactorInitNone<dim_outer, dim_tridiag> no_init_outer_subdiag;
    const FactorInitNone<dim_outer, dim_outer> no_init_outer_diag;
    fillWithOriginalMatrixValuesIncludingDiagonal(input.outer.diag,
                                                  m_data.outer.diag);
    ::ctldl::factorizePartialUpLooking(
        m_data.start.diag, m_data.start.next, input.start.outer,
        input.outer.subdiag[0], no_init_outer_diag, m_data.start.outer,
        m_data.outer.subdiag[0], m_data.outer.diag);
    for (std::size_t i = 0; i < num_repetitions; ++i) {
      ::ctldl::factorizePartialUpLooking(
          m_data.tridiag.diag[i], m_data.tridiag.subdiag[i],
          no_init_outer_subdiag, input.outer.subdiag[i + 1], no_init_outer_diag,
          m_data.outer.subdiag[i], m_data.outer.subdiag[i + 1],
          m_data.outer.diag);
    }
    ::ctldl::factorizePartialUpLooking(m_data.tridiag.diag[num_repetitions],
                                       m_data.link.prev, no_init_outer_subdiag,
                                       input.link.next, no_init_outer_diag,
                                       m_data.outer.subdiag[num_repetitions],
                                       m_data.link.next, m_data.outer.diag);
    ::ctldl::factorize(m_data.link.diag, FactorInitNone<dim_outer, dim_link>{},
                       no_init_outer_diag, m_data.link.next, m_data.outer.diag,
                       regularization, method_tag);
  }

  template <class Rhs>
  void solveInPlace(Rhs& rhs) const {
    forwardSolve(rhs);
    diagonalSolve(rhs);
    backwardSolve(rhs);
  }

  template <class Rhs>
  [[gnu::flatten]] void forwardSolve(Rhs& rhs) const {
    checkNumRepetitionsClassInvariant();

    const auto num_repetitions = std::size_t{numRepetitions()};
    solveForwardSubstitution(m_data.start.diag, rhs.start);
    solveForwardSubstitution(m_data.tridiag.diag[0], rhs.tridiag[0],
                             m_data.start.next, rhs.start);
    for (std::size_t i = 1; i <= num_repetitions; ++i) {
      solveForwardSubstitution(m_data.tridiag.diag[i], rhs.tridiag[i],
                               m_data.tridiag.subdiag[i - 1],
                               rhs.tridiag[i - 1]);
    }
    solveForwardSubstitution(m_data.link.diag, rhs.link, m_data.link.prev,
                             rhs.tridiag[num_repetitions]);
    solveForwardSubstitution(m_data.start.outer, rhs.start, rhs.outer);
    for (std::size_t i = 0; i <= num_repetitions; ++i) {
      solveForwardSubstitution(m_data.outer.subdiag[i], rhs.tridiag[i],
                               rhs.outer);
    }
    solveForwardSubstitution(m_data.link.next, rhs.link, rhs.outer);
    solveForwardSubstitution(m_data.outer.diag, rhs.outer);
  }

  template <class Rhs>
  void diagonalSolve(Rhs& rhs) const {
    checkNumRepetitionsClassInvariant();

    const auto num_repetitions = std::size_t{numRepetitions()};
    m_data.start.diag.diagonalSolve(rhs.start);
    for (std::size_t i = 0; i <= num_repetitions; ++i) {
      m_data.tridiag.diag[i].diagonalSolve(rhs.tridiag[i]);
    }
    m_data.link.diag.diagonalSolve(rhs.link);
    m_data.outer.diag.diagonalSolve(rhs.outer);
  }

  template <class Rhs>
  [[gnu::flatten]] void backwardSolve(Rhs& rhs) const {
    checkNumRepetitionsClassInvariant();

    const auto num_repetitions = std::size_t{numRepetitions()};
    // finish rhs_outer
    solveBackwardSubstitution(m_data.outer.diag, rhs.outer);
    // finish rhs_link
    solveBackwardSubstitution(m_data.link.next, rhs.outer, rhs.link);
    solveBackwardSubstitution(m_data.link.diag, rhs.link);
    // finish rhs[num_repetitions]
    solveBackwardSubstitution(m_data.outer.subdiag[num_repetitions],
                              rhs.outer, rhs.tridiag[num_repetitions]);
    solveBackwardSubstitution(m_data.link.prev, rhs.link,
                              rhs.tridiag[num_repetitions]);
    solveBackwardSubstitution(m_data.tridiag.diag[num_repetitions],
                              rhs.tridiag[num_repetitions]);
    // now everything else
    for (std::size_t i = num_repetitions; i > 0; --i) {
      solveBackwardSubstitution(m_data.outer.subdiag[i - 1], rhs.outer,
                                rhs.tridiag[i - 1]);
      solveBackwardSubstitution(m_data.tridiag.subdiag[i - 1], rhs.tridiag[i],
                                rhs.tridiag[i - 1]);
      solveBackwardSubstitution(m_data.tridiag.diag[i - 1], rhs.tridiag[i - 1]);
    }
    solveBackwardSubstitution(m_data.start.outer, rhs.outer, rhs.start);
    solveBackwardSubstitution(m_data.start.next, rhs.tridiag[0], rhs.start);
    solveBackwardSubstitution(m_data.start.diag, rhs.start);
  }

 private:
  MatrixTridiagonalArrowheadLinked<FactorStart, FactorTridiagonal, FactorLink,
                                   FactorOuter>
      m_data;
};

}  // namespace ctldl
