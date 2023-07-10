#pragma once

#include <ctldl/fileio/mtx_check.hpp>
#include <ctldl/fileio/mtx_foreach_entry.hpp>
#include <ctldl/fileio/mtx_read_header.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace ctldl {

template <class MatrixA, class MatrixB>
std::pair<std::vector<MatrixA>, std::vector<MatrixB>>
mtxFileReadRepeatingBlockTridiagonal(const char* filepath) {
  static_assert(MatrixA::sparsity.num_rows == MatrixA::sparsity.num_cols);
  static_assert(MatrixB::sparsity.num_rows == MatrixB::sparsity.num_cols);
  static_assert(MatrixA::sparsity.num_rows == MatrixB::sparsity.num_rows);
  constexpr auto dim = std::size_t{MatrixA::sparsity.num_rows};

  const std::size_t num_repetitions = [filepath] {
    const auto header = mtxReadHeader(filepath);
    mtxCheck(header.num_rows == header.num_cols, "Matrix must be square");
    mtxCheck(header.num_rows > 0, "Matrix dimension must be non-zero");
    mtxCheck(header.num_rows % dim == 0,
             "Matrix dimension must be multiple of expected repeating block "
             "dimension");
    return (header.num_rows / dim) - 1;
  }();

  std::vector<std::array<double, MatrixA::nnz>> values_A(num_repetitions + 1,
                                                         {0.0});
  std::vector<std::array<double, MatrixB::nnz>> values_B(num_repetitions,
                                                         {0.0});
  mtxForeachEntry(filepath, [&values_A, &values_B](const MtxEntry entry) {
    mtxCheck(entry.row_index >= entry.col_index,
             "Matrix entries must be in lower triangle");
    const std::size_t repetition_index = (entry.col_index / dim);
    const bool is_diagonal_block =
        (entry.row_index / dim == entry.col_index / dim);

    const auto block_row_index = entry.row_index % dim;
    const auto block_col_index = entry.col_index % dim;
    if (is_diagonal_block) {
      mtxCheck(MatrixA::sparsity.isNonZero(block_row_index, block_col_index),
               "Entry is not covered by compiled diagonal block sparsity");
      const auto entry_index =
          MatrixA::sparsity.entryIndex(block_row_index, block_col_index);
      values_A[repetition_index][entry_index] = entry.value;
    } else {
      mtxCheck(MatrixB::sparsity.isNonZero(block_row_index, block_col_index),
               "Entry is not covered by compiled subdiagonal block sparsity");
      const auto entry_index =
          MatrixB::sparsity.entryIndex(block_row_index, block_col_index);
      values_B[repetition_index][entry_index] = entry.value;
    }
  });

  std::vector<MatrixA> matrices_A(values_A.size(), MatrixA({0.0}));
  std::transform(values_A.cbegin(), values_A.cend(), matrices_A.begin(),
                 [](const auto& values) { return MatrixA{values}; });
  std::vector<MatrixB> matrices_B(values_B.size(), MatrixB({0.0}));
  std::transform(values_B.cbegin(), values_B.cend(), matrices_B.begin(),
                 [](const auto& values) { return MatrixB{values}; });

  return {matrices_A, matrices_B};
}

}  // namespace ctldl
