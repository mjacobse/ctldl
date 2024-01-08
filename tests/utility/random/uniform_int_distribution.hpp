#pragma once

#include <cassert>
#include <type_traits>

namespace ctldl {

template <class T>
class UniformIntDistribution {
 public:
  constexpr explicit UniformIntDistribution(const T a, const T b)
      : m_a(a), m_b(b) {
    assert(b >= a);
  }

  template <class UniformRandomNumberGenerator>
  constexpr T operator()(UniformRandomNumberGenerator& generator) {
    if (m_b == m_a) {
      return m_a;
    }
    const auto value = generator() - generator.min();
    const auto diff = static_cast<std::make_unsigned_t<T>>(m_b - m_a);
    // using modulo like this is bad since it is not uniform, but it is simpler
    // for now
    return static_cast<T>(value % diff) + m_a;
  }

 private:
  T m_a;
  T m_b;
};

}  // namespace ctldl
