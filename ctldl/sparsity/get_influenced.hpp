#pragma once

#include <array>
#include <cstddef>

namespace ctldl {

template <class SparsitySource>
constexpr std::size_t getNumInfluenced(const std::size_t j,
                                       const std::size_t k_end) {
  std::size_t num = 0;
  for (std::size_t k = 0; k < k_end; ++k) {
    if (SparsitySource::isNonZero(k, j)) {
      num += 1;
    }
  }
  return num;
};

struct Influenced {
  std::size_t entry_index_source;
  std::size_t entry_index_target;
};

template <std::size_t i, std::size_t j, class SparsitySource,
          class SparsityTarget, class Limiter>
constexpr auto getInfluencedList(Limiter limiter) {
  constexpr auto k_end = limiter(i, j);
  constexpr auto num_influenced = getNumInfluenced<SparsitySource>(j, k_end);
  std::array<Influenced, num_influenced> influenced_list;
  std::size_t influenced_index = 0;
  for (std::size_t k = 0; k < k_end; ++k) {
    if (SparsitySource::isNonZero(k, j)) {
      influenced_list[influenced_index] = Influenced{
          SparsitySource::entryIndex(k, j), SparsityTarget::entryIndex(i, k)};
      influenced_index += 1;
    }
  }
  return influenced_list;
}

}  // namespace ctldl
