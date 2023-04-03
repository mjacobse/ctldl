#pragma once

#include <ctldl/sparsity/get_entries.hpp>
#include <ctldl/sparsity/get_filled_in_is_nonzero_info.hpp>
#include <ctldl/sparsity/get_is_nonzero_info.hpp>

#include <cassert>

namespace ctldl {

template <class Sparsity>
struct FilledInSparsity {
  static constexpr auto num_rows = Sparsity::num_rows;
  static constexpr auto num_cols = Sparsity::num_cols;
  static constexpr auto is_nonzero =
      getFilledInIsNonzeroInfo(getIsNonzeroInfo<Sparsity>());
  static constexpr auto entries = getEntries([] { return is_nonzero; });
};

}  // namespace ctldl
