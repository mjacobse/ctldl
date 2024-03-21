#pragma once

namespace ctldl {

template <class Start, class Tridiag, class Link, class Outer>
struct VectorTridiagonalArrowheadLinked {
  Start start;
  Tridiag tridiag;
  Link link;
  Outer outer;
};

// needed for clang < 17 which does not do the CTAD for aggregates otherwise
template <class Start, class Tridiag, class Link, class Outer>
VectorTridiagonalArrowheadLinked(Start, Tridiag, Link, Outer)
    -> VectorTridiagonalArrowheadLinked<Start, Tridiag, Link, Outer>;

}  // namespace ctldl
