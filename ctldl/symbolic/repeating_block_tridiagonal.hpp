#pragma once

namespace ctldl {

template <class A, class B>
struct RepeatingBlockTridiagonal {
  using Diagonal = A;
  using Subdiagonal = B;
  Diagonal diagonal;
  Subdiagonal subdiagonal;
};

template <class A, class B>
RepeatingBlockTridiagonal(A, B) -> RepeatingBlockTridiagonal<A, B>;

}  // namespace ctldl
