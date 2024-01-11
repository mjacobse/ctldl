#pragma once

namespace ctldl {

template <class Diag, class Subdiag>
struct MatrixTridiagonal {
  Diag diag;
  Subdiag subdiag;
};

// needed for clang < 17 which does not do the CTAD for aggregates otherwise
template <class Diag, class Subdiag>
MatrixTridiagonal(Diag, Subdiag) -> MatrixTridiagonal<Diag, Subdiag>;

}  // namespace ctldl
