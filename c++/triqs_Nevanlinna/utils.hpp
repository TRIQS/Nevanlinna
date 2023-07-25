#ifndef TRIQS_NEVANLINNA_UTILS_HPP
#define TRIQS_NEVANLINNA_UTILS_HPP

#include <cstdlib>
#include <string>

namespace triqs_Nevanlinna {
  /**
   * Extract and return numerical value of environment variable.
   * If the environment variable is not set return the defult value.
   *
   * @param env_var - environment variable name
   * @param def - defualt value
   * @return
   */
  inline int get_env_int(const std::string& env_var, int def) {
    if(const char * val = std::getenv(env_var.c_str())) {
      return std::atoi(val);
    }
    return def;
  }
}

#endif //TRIQS_NEVANLINNA_UTILS_HPP
