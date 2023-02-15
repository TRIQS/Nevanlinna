#ifndef NEVANLINNA_NEVANLINNA_ERROR_HPP
#define NEVANLINNA_NEVANLINNA_ERROR_HPP

#include <stdexcept>

namespace triqs_Nevanlinna {
  class Nevanlinna_error : public std::runtime_error {
    public:
    explicit Nevanlinna_error(const std::string &msg) : std::runtime_error(msg) {}
  };

  class Nevanlinna_uninitialized_error : public Nevanlinna_error {
    public:
    explicit Nevanlinna_uninitialized_error(const std::string &msg) : Nevanlinna_error(msg) {}
  };

  class Nevanlinna_negative_grid_error : public Nevanlinna_error {
    public:
    explicit Nevanlinna_negative_grid_error(const std::string &msg) : Nevanlinna_error(msg) {}
  };

  class unimplemented_kernel_error : public Nevanlinna_error {
    public:
    explicit unimplemented_kernel_error(const std::string &msg) : Nevanlinna_error(msg) {}
  };

} // namespace triqs_Nevanlinna
#endif //NEVANLINNA_NEVANLINNA_ERROR_HPP
