#pragma once

#include <array>
#include <cstddef>

namespace ctldl {

template <std::size_t num_rows_, std::size_t num_cols_>
class IsNonzeroInfo {
 public:
  static constexpr auto num_rows = num_rows_;
  static constexpr auto num_cols = num_cols_;

  constexpr std::array<bool, num_cols>& operator[](const std::size_t i) {
    return m_data[i];
  }

  constexpr const std::array<bool, num_cols>& operator[](
      const std::size_t i) const {
    return m_data[i];
  }

  constexpr std::size_t nnz() const {
    std::size_t count = 0;
    for (std::size_t i = 0; i < num_rows; ++i) {
      for (std::size_t j = 0; j < num_cols; ++j) {
        count += m_data[i][j];
      }
    }
    return count;
  }

 private:
  std::array<std::array<bool, num_cols>, num_rows> m_data{};
};

}  // namespace ctldl
