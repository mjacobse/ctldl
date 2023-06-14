#include <ctldl/fileio/mtx_read_header.hpp>

#include <ctldl/fileio/mtx_check.hpp>

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <string>

namespace ctldl {

MtxDimensions mtxReadHeader(const char* filepath) {
  std::ifstream file(filepath);
  return mtxReadHeader(file);
}

MtxDimensions mtxReadHeader(std::ifstream& file) {
  std::string line;
  while (std::getline(file, line) && line[0] == '%') {
    // skip comment lines
  }

  std::size_t num_rows;
  std::size_t num_cols;
  std::size_t nnz;
  const auto num_read =
      std::sscanf(line.c_str(), " %zu %zu %zu", &num_rows, &num_cols, &nnz);
  mtxCheck(num_read == 3, "Error reading matrix dimensions");
  return {num_rows, num_cols, nnz};
}

}  // namespace ctldl
