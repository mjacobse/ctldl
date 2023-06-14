#include <ctldl/fileio/mtx_check.hpp>

#include <stdexcept>

namespace ctldl {

class MtxReadException : public std::runtime_error {
 public:
  explicit MtxReadException(const char* error_message)
      : std::runtime_error(error_message) {}
};

void mtxCheck(const bool condition, const char* message) {
  if (!condition) {
    throw MtxReadException(message);
  }
}

}  // namespace ctldl
