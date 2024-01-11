#pragma once

namespace ctldl {

template <class Tridiag, class Link, class Outer>
struct VectorTridiagonalArrowheadLinked {
  Tridiag tridiag;
  Link link;
  Outer outer;
};

// needed for clang < 17 which does not do the CTAD for aggregates otherwise
template <class Tridiag, class Link, class Outer>
VectorTridiagonalArrowheadLinked(Tridiag, Link, Outer)
    -> VectorTridiagonalArrowheadLinked<Tridiag, Link, Outer>;

}  // namespace ctldl
