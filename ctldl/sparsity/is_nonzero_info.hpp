#pragma once

#include <array>

namespace ctldl {

template <int num_rows_, int num_cols_>
class IsNonzeroInfo {
 public:
  static constexpr auto num_rows = num_rows_;
  static constexpr auto num_cols = num_cols_;

  constexpr std::array<bool, num_cols>& operator[](const int i) {
    return m_data[i];
  }

  constexpr const std::array<bool, num_cols>& operator[](const int i) const {
    return m_data[i];
  }

  constexpr int nnz() const {
    int count = 0;
    for (int i = 0; i < num_rows; ++i) {
      for (int j = 0; j < num_cols; ++j) {
        count += m_data[i][j];
      }
    }
    return count;
  }

 private:
  std::array<std::array<bool, num_cols>, num_rows> m_data{false};
};

}  // namespace ctldl
