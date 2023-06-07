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
  } else if constexpr (method == FactorizeMethod::UpLooking) {
    return FactorizeMethodUpLooking{};
  } else {
    // Basically static_assert(false), but need to make it dependent on template
    // parameter, otherwise it always fires.
    static_assert(method != method, "Invalid FactorizeMethod.");
  }
}

}  // namespace ctldl
