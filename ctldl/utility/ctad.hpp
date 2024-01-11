#pragma once

#include <utility>

namespace ctldl {

// Gives the type for class template T that is deduced through class template
// argument deduction for arguments of type U.
// Only works for class templates that take exclusively non-type template
// parameters for now (not easy to write a general solution before P2989).
template <template <auto...> class T, class U>
struct ctad {
  using type = decltype(T{std::declval<U>()});
};

template <template <auto...> class T, class U>
using ctad_t = typename ctad<T, U>::type;

}  // namespace ctldl
