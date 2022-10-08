//
// Created by iskakoff on 10/7/22.
//

#ifndef NEVANLINNA_NEVANLINNA_ERROR_HPP
#define NEVANLINNA_NEVANLINNA_ERROR_HPP

#include <stdexcept>

namespace nevanlinna {
  class nevanlinna_error : public std::runtime_error {
    public:
    explicit
    nevanlinna_error(const std::string &msg) : std::runtime_error(msg) {}

  };

}
#endif //NEVANLINNA_NEVANLINNA_ERROR_HPP
