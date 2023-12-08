#pragma once

namespace ctldl {

template <class Result, class Generator>
struct CompileTimeRandomResult {
  Result result;
  Generator generator;
};

template <class Distribution, class Generator>
constexpr auto compileTimeGenerate(const Distribution& distribution,
                                   Generator generator) {
  auto result = distribution(generator);
  return CompileTimeRandomResult<decltype(result), Generator>{result,
                                                              generator};
}

}  // namespace ctldl
