#pragma once

namespace ctldl {

template <class Diagonal, class Subdiagonal, class Outer>
struct RepeatingBlockTridiagonalArrowhead {
  Diagonal diag;
  Subdiagonal subdiag;
  Outer outer;
};

template <class Diagonal, class Subdiagonal, class Outer>
RepeatingBlockTridiagonalArrowhead(Diagonal, Subdiagonal, Outer)
    -> RepeatingBlockTridiagonalArrowhead<Diagonal, Subdiagonal, Outer>;

}  // namespace ctldl
