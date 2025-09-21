#pragma once

namespace ctldl {

template <class Diag, class Subdiag>
struct MatrixTridiagonal {
  Diag diag;
  Subdiag subdiag;
};

}  // namespace ctldl
