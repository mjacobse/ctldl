#pragma once

namespace ctldl {

enum class FactorizeMethod {
    EntryWise,
    UpLooking,
};

struct FactorizeMethodEntryWise{};
struct FactorizeMethodUpLooking{};

template <FactorizeMethod method>
auto getFactorizeMethodTag() {
  if constexpr (method == FactorizeMethod::EntryWise) {
    return FactorizeMethodEntryWise{};
  } else {
    static_assert(method == FactorizeMethod::UpLooking,
                  "Invalid FactorizeMethod.");
    return FactorizeMethodUpLooking{};
  }
}

}  // namespace ctldl
