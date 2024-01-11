#pragma once

#include <cassert>
#include <cstddef>

namespace ctldl {

template <class T>
class DummyRange {
 public:
  explicit DummyRange(const std::size_t n) : m_size(n) {};

  constexpr T operator[](const std::size_t i) const {
    assert(i < m_size);
    static_cast<void>(i);
    return T{};
  }

 private:
  std::size_t m_size;
};

}  // namespace ctldl
