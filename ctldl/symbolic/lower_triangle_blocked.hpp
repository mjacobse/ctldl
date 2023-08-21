#pragma once

namespace ctldl {

template <class T11, class T21, class T22>
struct LowerTriangleBlocked {
  using Block11 = T11;
  using Block21 = T21;
  using Block22 = T22;
  Block11 block11;
  Block21 block21;
  Block22 block22;
};

template <class T11, class T21, class T22>
LowerTriangleBlocked(T11, T21, T22) -> LowerTriangleBlocked<T11, T21, T22>;

}  // namespace ctldl
