#include <ctldl/fileio/mtx_foreach_entry.hpp>

#include <ctldl/fileio/mtx_check.hpp>
#include <ctldl/fileio/mtx_read_header.hpp>

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>

namespace ctldl {

void mtxForeachEntry(const char* filepath, std::function<void(MtxEntry)> func) {
  constexpr std::size_t row_index_base = 1;
  constexpr std::size_t col_index_base = 1;

  std::ifstream file(filepath);
  const auto dimensions = mtxReadHeader(file);

  std::size_t num_entries = 0;
  std::string line;
  while (std::getline(file, line)) {
    std::size_t row_index;
    std::size_t col_index;
    double value;
    const auto num_read = std::sscanf(line.c_str(), " %zu %zu %lf",
                                      &row_index, &col_index, &value);
    mtxCheck(num_read == 3, "Error reading matrix entry");
    mtxCheck(row_index > 0, "Row index must be positive");
    mtxCheck(col_index > 0, "Column index must be positive");
    mtxCheck(row_index <= dimensions.num_rows, "Row index is out of range");
    mtxCheck(col_index <= dimensions.num_cols, "Column index is out of range");
    func(MtxEntry{row_index - row_index_base, col_index - col_index_base,
                  value});
    num_entries += 1;
  }
  mtxCheck(
      num_entries == dimensions.nnz,
      "Number of nonzeros in header does not equal number of read entries");
}

}  // namespace ctldl
