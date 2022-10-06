#pragma once
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
//#include <triqs/arrays.hpp>
#include <nda/nda.hpp>
#include <h5/h5.hpp>
#include <gmpxx.h>

#include "nevanlinna_parameters_t.hpp"

namespace nevanlinna {

  /**
   * Nevanlinna analytical continuation solver
   *
   * @note A Useful note
   * @include nevanlinna/nevanlinna.hpp
   */
  //template<typename T>
  class solver {

    public:
//    solver() = default;
    ~solver() = default;

    CPP2PY_ARG_AS_DICT
    solver(nevanlinna_parameters_t const & p);

    // Copy/Move construction
    solver(solver const &) = default;
    solver(solver &&)      = default;

    /// Copy/Move assignment
    solver &operator=(solver const &) = default;
    solver &operator=(solver &&) = default;

    void solve(const nda::array<double, 2> & input);

    nda::array<double, 1> evaluate(const nda::array<double, 1> & grid);
  };

//  typedef solver<mpf_class> SolverMpf;
} // namespace nevanlinna
