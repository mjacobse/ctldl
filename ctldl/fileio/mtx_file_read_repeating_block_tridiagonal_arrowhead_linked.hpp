#pragma once

#include <ctldl/fileio/mtx_check.hpp>
#include <ctldl/fileio/mtx_foreach_entry.hpp>
#include <ctldl/fileio/mtx_read_header.hpp>
#include <ctldl/matrix/matrix_link.hpp>
#include <ctldl/matrix/matrix_outer.hpp>
#include <ctldl/matrix/matrix_start.hpp>
#include <ctldl/matrix/matrix_tridiagonal.hpp>
#include <ctldl/matrix/matrix_tridiagonal_arrowhead_linked.hpp>
#include <ctldl/sparsity/is_square.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace ctldl {

template <class Matrix>
void mtxAddEntry(const MtxEntry entry, Matrix& matrix) {
  mtxCheck(Matrix::sparsity.isNonZero(entry.row_index, entry.col_index),
           "Entry is not covered by compiled sparsity");
  const auto entry_index =
      Matrix::sparsity.entryIndex(entry.row_index, entry.col_index);
  matrix.values[entry_index] = entry.value;
}

template <auto sparsity, class Matrix>
void mtxAddEntryStart(const MtxEntry entry, Matrix& matrix) {
  if (entry.row_index < sparsity.dim) {
    mtxAddEntry(entry, matrix.diag);
  } else {
    mtxAddEntry({entry.row_index - sparsity.dim, entry.col_index, entry.value},
                matrix.next);
  }
}

template <auto sparsity, class Matrix>
void mtxAddEntryTridiag(const MtxEntry entry, Matrix& matrix) {
  const std::size_t repetition_index = (entry.col_index / sparsity.dim);
  const bool is_diagonal_block =
      (entry.row_index / sparsity.dim == entry.col_index / sparsity.dim);

  const auto block_entry =
      MtxEntry{entry.row_index % sparsity.dim, entry.col_index % sparsity.dim,
               entry.value};
  if (is_diagonal_block) {
    mtxAddEntry(block_entry, matrix.diag[repetition_index]);
  } else {
    mtxAddEntry(block_entry, matrix.subdiag[repetition_index]);
  }
}

template <auto sparsity, class Matrix>
void mtxAddEntryOuterSubdiag(const MtxEntry entry, Matrix& matrix) {
  const std::size_t repetition_index = (entry.col_index / sparsity.dim_inner);

  const auto block_entry =
      MtxEntry{entry.row_index % sparsity.dim_inner,
               entry.col_index % sparsity.dim_inner, entry.value};
  mtxAddEntry(block_entry, matrix.subdiag[repetition_index]);
}

template <auto sparsity, class Matrix>
void mtxAddEntryLink(const MtxEntry entry, Matrix& matrix) {
  if (entry.col_index < sparsity.dim_prev) {
    mtxAddEntry(entry, matrix.prev);
    return;
  }

  if (entry.row_index >= sparsity.dim) {
    mtxAddEntry({entry.row_index - sparsity.dim,
                 entry.col_index - sparsity.dim_prev, entry.value},
                matrix.next);
    return;
  }

  mtxAddEntry(
      {entry.row_index, entry.col_index - sparsity.dim_prev, entry.value},
      matrix.diag);
}

template <auto sparsity_in>
struct MtxInputMatrix {
  static constexpr auto sparsity = SparsityCSR(sparsity_in);
  std::array<double, sparsity.nnz> values = {};
  double valueAt(const std::size_t entry_index) const {
    return values[entry_index];
  }
};

template <auto sparsity>
auto mtxFileReadRepeatingBlockTridiagonalArrowheadLinked(const char* filepath) {
  const std::size_t num_repetitions = [filepath] {
    const auto header = mtxReadHeader(filepath);
    mtxCheck(header.num_rows == header.num_cols, "Matrix must be square");
    const auto num_rows_nonrepeating =
        sparsity.dim_start + sparsity.dim_link + sparsity.dim_outer;
    mtxCheck(header.num_rows > num_rows_nonrepeating,
             "Matrix dimension must be greater than sum of dimensions of "
             "non-repeating blocks.");
    const auto num_rows_tridiag = header.num_rows - num_rows_nonrepeating;
    mtxCheck(num_rows_tridiag % sparsity.dim_tridiag == 0,
             "Dimension of repeating part of matrix must be multiple of "
             "expected repeating block dimension");
    return (num_rows_tridiag / sparsity.dim_tridiag) - 1;
  }();

  MatrixStart start = {
      MtxInputMatrix<sparsity.start.diag>{},
      MtxInputMatrix<sparsity.start.next>{},
      MtxInputMatrix<sparsity.start.outer>{},
  };
  MatrixTridiagonal tridiag = {
      std::vector<MtxInputMatrix<sparsity.tridiag.diag>>(num_repetitions + 1),
      std::vector<MtxInputMatrix<sparsity.tridiag.subdiag>>(num_repetitions),
  };
  MatrixLink link = {
      MtxInputMatrix<sparsity.link.prev>{},
      MtxInputMatrix<sparsity.link.diag>{},
      MtxInputMatrix<sparsity.link.next>{},
  };
  MatrixOuter outer = {
      std::vector<MtxInputMatrix<sparsity.outer.subdiag>>(num_repetitions + 1),
      MtxInputMatrix<sparsity.outer.diag>{},
  };

  constexpr auto tridiag_begin = sparsity.dim_start;
  const auto tridiag_end =
      tridiag_begin + ((num_repetitions + 1) * sparsity.dim_tridiag);
  const auto outer_begin = tridiag_end + sparsity.dim_link;

  mtxForeachEntry(filepath, [tridiag_end, outer_begin, &start, &tridiag, &link,
                             &outer](const MtxEntry entry) {
    mtxCheck(entry.row_index >= entry.col_index,
             "Matrix entries must be in lower triangle");
    if (entry.col_index < tridiag_begin) {
      if (entry.row_index < outer_begin) {
        mtxAddEntryStart<sparsity.start>(entry, start);
      } else {
        mtxAddEntry(
            {entry.row_index - outer_begin, entry.col_index, entry.value},
            start.outer);
      }
      return;
    }

    if (entry.row_index < tridiag_end) {
      mtxAddEntryTridiag<sparsity.tridiag>(
          {entry.row_index - tridiag_begin, entry.col_index - tridiag_begin,
           entry.value},
          tridiag);
      return;
    }

    if (entry.row_index >= outer_begin && entry.col_index >= outer_begin) {
      mtxAddEntry({entry.row_index - outer_begin, entry.col_index - outer_begin,
                   entry.value},
                  outer.diag);
      return;
    }

    if (entry.row_index >= outer_begin && entry.col_index < tridiag_end) {
      mtxAddEntryOuterSubdiag<sparsity.outer>(
          {entry.row_index - outer_begin, entry.col_index - tridiag_begin,
           entry.value},
          outer);
      return;
    }

    const auto tridiag_last_block_begin = tridiag_end - sparsity.dim_tridiag;
    mtxCheck(entry.col_index >= tridiag_last_block_begin,
             "Entry is not covered by compiled sparsity");
    mtxAddEntryLink<sparsity.link>(
        {entry.row_index - tridiag_end,
         entry.col_index - tridiag_last_block_begin, entry.value},
        link);
  });

  return MatrixTridiagonalArrowheadLinked{start, tridiag, link, outer};
}

}  // namespace ctldl
