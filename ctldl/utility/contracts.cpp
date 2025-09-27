#include <ctldl/utility/contracts.hpp>

#include <format>
#include <stdexcept>
#include <string_view>

namespace ctldl {

class ContractViolationException : public std::logic_error {
 public:
  explicit ContractViolationException(
      const std::string_view kind, const std::source_location source_location =
                                       std::source_location::current())
      : std::logic_error(
            std::format("Contract {} violation at {}:{} in function {}", kind,
                        source_location.file_name(), source_location.line(),
                        source_location.function_name())) {}
};

void contractHandlePreconditionViolation(
    const std::source_location source_location) {
  throw ContractViolationException("precondition", source_location);
}

void contractHandlePostconditionViolation(
    const std::source_location source_location) {
  throw ContractViolationException("postcondition", source_location);
}

void contractHandleAssertionViolation(
    const std::source_location source_location) {
  throw ContractViolationException("assertion", source_location);
}

}  // namespace ctldl
