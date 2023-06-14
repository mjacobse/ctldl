#pragma once

#include <cstddef>
#include <functional>

namespace ctldl {

struct MtxEntry {
  std::size_t row_index;
  std::size_t col_index;
  double value;
};

void mtxForeachEntry(const char* filepath, std::function<void(MtxEntry)> func);

}  // namespace ctldl
