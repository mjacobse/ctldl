#pragma once

#include <cstddef>

namespace ctldl {

class BernoulliDistribution {
 public:
  constexpr explicit BernoulliDistribution(const double probability)
      : m_probability(probability) {}

  template <class UniformRandomNumberGenerator>
  constexpr bool operator()(UniformRandomNumberGenerator& generator) {
    // this is bad and should do something like generate_canonical() <
    // m_probability instead, but that would have to be reimplemented to be
    // usable in a constexpr context
    return (generator() - generator.min()) <
           static_cast<std::size_t>(
               m_probability *
               static_cast<double>(generator.max() - generator.min()));
  }

 private:
  double m_probability;
};

}  // namespace ctldl
