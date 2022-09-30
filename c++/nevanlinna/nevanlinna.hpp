#pragma once
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <h5/h5.hpp>

namespace nevanlinna {

  /**
   * Nevanlinna analytical continuation solver
   *
   * @note A Useful note
   * @include nevanlinna/nevanlinna.hpp
   */
  template<typename T>
  class solver {

    public:
    solver() = default;
    ~solver() = default;

    // Copy/Move construction
    solver(solver const &) = default;
    solver(solver &&)      = default;

    /// Copy/Move assignment
    solver &operator=(solver const &) = default;
    solver &operator=(solver &&) = default;
  };

} // namespace nevanlinna
