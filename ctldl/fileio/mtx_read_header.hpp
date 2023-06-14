#pragma once

#include <cstddef>
#include <iosfwd>

namespace ctldl {

struct MtxDimensions {
  std::size_t num_rows;
  std::size_t num_cols;
  std::size_t nnz;
};

MtxDimensions mtxReadHeader(const char* filepath);
MtxDimensions mtxReadHeader(std::ifstream& file);

}  // namespace ctldl
