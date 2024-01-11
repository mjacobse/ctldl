#pragma once

namespace ctldl {

template <class T11, class T21, class T22, class T31, class T32, class T33>
struct LowerTriangleBlocked3x3 {
  using Block11 = T11;
  using Block21 = T21;
  using Block22 = T22;
  using Block31 = T31;
  using Block32 = T32;
  using Block33 = T33;
  Block11 block11;
  Block21 block21;
  Block22 block22;
  Block31 block31;
  Block32 block32;
  Block33 block33;
};

// needed for clang < 17 which does not do the CTAD for aggregates otherwise
template <class T11, class T21, class T22, class T31, class T32, class T33>
LowerTriangleBlocked3x3(T11, T21, T22, T31, T32, T33)
    -> LowerTriangleBlocked3x3<T11, T21, T22, T31, T32, T33>;

template <class T11, class T21, class T22>
struct LowerTriangleBlocked {
  using Block11 = T11;
  using Block21 = T21;
  using Block22 = T22;
  Block11 block11;
  Block21 block21;
  Block22 block22;
};

// needed for clang < 17 which does not do the CTAD for aggregates otherwise
template <class T11, class T21, class T22>
LowerTriangleBlocked(T11, T21, T22) -> LowerTriangleBlocked<T11, T21, T22>;

}  // namespace ctldl
