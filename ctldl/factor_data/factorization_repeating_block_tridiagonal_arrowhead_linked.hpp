#pragma once

#include <ctldl/factor_data/factorization_already_permuted.hpp>
#include <ctldl/factor_data/factorization_subdiagonal_block.hpp>
#include <ctldl/factor_data/sparsity_to_factorize_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/factorize/factorize_entry_wise.hpp>
#include <ctldl/factorize/factorize_method.hpp>
#include <ctldl/factorize/factorize_up_looking.hpp>
#include <ctldl/permutation/permutation.hpp>
#include <ctldl/solve/solve_backward_substitution.hpp>
#include <ctldl/solve/solve_forward_substitution.hpp>
#include <ctldl/sparsity/sparsity_lower_triangle.hpp>
#include <ctldl/sparsity/sparsity_permuted.hpp>
#include <ctldl/symbolic/filled_in_sparsity_blocked.hpp>
#include <ctldl/symbolic/filled_in_sparsity_repeating.hpp>

#include <cstddef>
#include <memory>
#include <span>


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
  static constexpr auto dim_tridiag = std::size_t{sparsity.dim_tridiag};
  static constexpr auto dim_link = std::size_t{sparsity.dim_link};
  static constexpr auto dim_outer = std::size_t{sparsity.dim_outer};
  static constexpr auto permutation_tridiag = sparsity.tridiag.permutation;
  static constexpr auto permutation_link = sparsity.link.permutation;
  static constexpr auto permutation_outer = sparsity.outer.permutation;

  static constexpr auto sparsity_factor = getFilledInSparsityRepeatingArrowhead<
      sparsity.tridiag.diag, sparsity.tridiag.subdiag, sparsity.outer.subdiag,
      permutation_tridiag, permutation_outer>();

  static constexpr auto helper = [] {
    constexpr auto sparsity_link_tridiag_permuted_cols = getSparsityPermuted(
        sparsity.link.prev, Permutation<dim_link>{}, permutation_tridiag);
    constexpr auto sparsity_link_outer_permuted_rows = getSparsityPermuted(
        sparsity.link.next, permutation_outer, Permutation<dim_link>{});
    constexpr auto sparsity_outer_diag_permuted =
        getSparsityLowerTriangle<sparsity.outer.diag>(permutation_outer);
    return getFilledInSparsityBlocked3x3<
        sparsity_factor.diag, sparsity_link_tridiag_permuted_cols,
        sparsity.link.diag, sparsity_factor.outer,
        sparsity_link_outer_permuted_rows, sparsity_outer_diag_permuted,
        Permutation<dim_tridiag>{}, permutation_link,
        Permutation<dim_outer>{}>();
  }();
  static constexpr auto sparsity_factor_link_tridiag = helper.block21;
  static constexpr auto sparsity_factor_link_diag = helper.block22;
  static constexpr auto sparsity_factor_link_outer = helper.block32;
  static constexpr auto sparsity_factor_outer_diag = helper.block33;

 public:
  using Value = Value_;
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

  explicit FactorizationRepeatingBlockTridiagonalArrowheadLinked(
      const std::size_t num_repetitions)
      : m_num_repetitions(num_repetitions),
        m_tridiag_diag(new FactorTridiagDiag[num_repetitions + 1]),
        m_tridiag_subdiag(new FactorTridiagSubdiag[num_repetitions]),
        m_outer_subdiag(new FactorOuterSubdiag[num_repetitions + 1]) {}

  auto numRepetitions() const noexcept { return m_num_repetitions; }
  std::span<const FactorTridiagDiag> blocksA() const noexcept {
    return {m_tridiag_diag.get(), m_num_repetitions + 1};
  }
  std::span<const FactorTridiagSubdiag> blocksB() const noexcept {
    return {m_tridiag_subdiag.get(), m_num_repetitions};
  }

  template <class FactorizeMethodTag = FactorizeMethodUpLooking,
            class MatrixBlocksInputValues>
  void factorize(const MatrixBlocksInputValues& input,
                 const FactorizeMethodTag method_tag = {}) {
    // tridiagonal
    m_tridiag_diag[0].factorize(input.tridiag.diag[0], method_tag);
    for (std::size_t i = 0; i < m_num_repetitions; ++i) {
      ::ctldl::factorize(m_tridiag_diag[i], input.tridiag.subdiag[i],
                         input.tridiag.diag[i + 1], m_tridiag_subdiag[i],
                         m_tridiag_diag[i + 1], method_tag);
    }
    ::ctldl::factorize(m_tridiag_diag[m_num_repetitions], input.link.prev,
                       input.link.diag, m_link_tridiag, m_link_diag,
                       method_tag);
    // outer
    const FactorInitNone<dim_outer, dim_tridiag> no_init_outer_subdiag;
    const FactorInitNone<dim_outer, dim_outer> no_init_outer_diag;
    fillWithOriginalMatrixValuesIncludingDiagonal(input.outer.diag,
                                                  m_outer_diag);
    fillWithOriginalMatrixValues(input.outer.subdiag[0], m_outer_subdiag[0]);
    for (std::size_t i = 0; i < m_num_repetitions; ++i) {
      ::ctldl::factorizePartialUpLooking(
          m_tridiag_diag[i], m_tridiag_subdiag[i], no_init_outer_subdiag,
          input.outer.subdiag[i + 1], no_init_outer_diag, m_outer_subdiag[i],
          m_outer_subdiag[i + 1], m_outer_diag);
    }
    ::ctldl::factorizePartialUpLooking(
        m_tridiag_diag[m_num_repetitions], m_link_tridiag,
        no_init_outer_subdiag, input.link.next, no_init_outer_diag,
        m_outer_subdiag[m_num_repetitions], m_link_outer, m_outer_diag);
    ::ctldl::factorize(m_link_diag, FactorInitNone<dim_outer, dim_link>{},
                       no_init_outer_diag, m_link_outer, m_outer_diag,
                       method_tag);
  }

  template <class Rhs>
  void solveInPlace(Rhs& rhs) const {
    forwardSolve(rhs);
    diagonalSolve(rhs);
    backwardSolve(rhs);
  }

  template <class Rhs>
  [[gnu::flatten]] void forwardSolve(Rhs& rhs) const {
    solveForwardSubstitution(m_tridiag_diag[0], rhs.tridiag[0]);
    for (std::size_t i = 1; i <= m_num_repetitions; ++i) {
      solveForwardSubstitution(m_tridiag_diag[i], rhs.tridiag[i],
                               m_tridiag_subdiag[i - 1], rhs.tridiag[i - 1]);
    }
    solveForwardSubstitution(m_link_diag, rhs.link, m_link_tridiag,
                             rhs.tridiag[m_num_repetitions]);
    for (std::size_t i = 0; i <= m_num_repetitions; ++i) {
      solveForwardSubstitution(m_outer_subdiag[i], rhs.tridiag[i], rhs.outer);
    }
    solveForwardSubstitution(m_link_outer, rhs.link, rhs.outer);
    solveForwardSubstitution(m_outer_diag, rhs.outer);
  }

  template <class Rhs>
  void diagonalSolve(Rhs& rhs) const {
    for (std::size_t i = 0; i <= m_num_repetitions; ++i) {
      m_tridiag_diag[i].diagonalSolve(rhs.tridiag[i]);
    }
    m_link_diag.diagonalSolve(rhs.link);
    m_outer_diag.diagonalSolve(rhs.outer);
  }

  template <class Rhs>
  [[gnu::flatten]] void backwardSolve(Rhs& rhs) const {
    // finish rhs_outer
    solveBackwardSubstitution(m_outer_diag, rhs.outer);
    // finish rhs_link
    solveBackwardSubstitution(m_link_outer, rhs.outer, rhs.link);
    solveBackwardSubstitution(m_link_diag, rhs.link);
    // finish rhs[m_num_repetitions]
    solveBackwardSubstitution(m_outer_subdiag[m_num_repetitions], rhs.outer,
                              rhs.tridiag[m_num_repetitions]);
    solveBackwardSubstitution(m_link_tridiag, rhs.link,
                              rhs.tridiag[m_num_repetitions]);
    solveBackwardSubstitution(m_tridiag_diag[m_num_repetitions],
                              rhs.tridiag[m_num_repetitions]);
    // now everything else
    for (std::size_t i = m_num_repetitions; i > 0; --i) {
      solveBackwardSubstitution(m_outer_subdiag[i - 1], rhs.outer,
                                rhs.tridiag[i - 1]);
      solveBackwardSubstitution(m_tridiag_subdiag[i - 1], rhs.tridiag[i],
                                rhs.tridiag[i - 1]);
      solveBackwardSubstitution(m_tridiag_diag[i - 1], rhs.tridiag[i - 1]);
    }
  }

 private:
  std::size_t m_num_repetitions;
  std::unique_ptr<FactorTridiagDiag[]> m_tridiag_diag;
  std::unique_ptr<FactorTridiagSubdiag[]> m_tridiag_subdiag;
  FactorLinkTridiag m_link_tridiag;
  FactorLinkDiag m_link_diag;
  FactorLinkOuter m_link_outer;
  std::unique_ptr<FactorOuterSubdiag[]> m_outer_subdiag;
  FactorOuterDiag m_outer_diag;
};

}  // namespace ctldl
