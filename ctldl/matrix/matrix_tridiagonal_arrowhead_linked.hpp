#pragma once

namespace ctldl {

template <class Start, class Tridiag, class Link, class Outer>
struct MatrixTridiagonalArrowheadLinked {
  Start start;
  Tridiag tridiag;
  Link link;
  Outer outer;
};

}  // namespace ctldl
