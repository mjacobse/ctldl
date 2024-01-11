#pragma once

namespace ctldl {

template <class Tridiag, class Link, class Outer>
struct MatrixTridiagonalArrowheadLinked {
  Tridiag tridiag;
  Link link;
  Outer outer;
};

// needed for clang < 17 which does not do the CTAD for aggregates otherwise
template <class Tridiag, class Link, class Outer>
MatrixTridiagonalArrowheadLinked(Tridiag, Link, Outer)
    -> MatrixTridiagonalArrowheadLinked<Tridiag, Link, Outer>;

}  // namespace ctldl
