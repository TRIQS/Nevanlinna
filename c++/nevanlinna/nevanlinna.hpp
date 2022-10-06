#pragma once
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <nda/nda.hpp>
#include <h5/h5.hpp>
#include <gmpxx.h>

#include "Eigen/Core"

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
    using complex_t = std::complex<mpf_class>;
    using matrix_t = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;

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

    void solve(const nda::array<std::complex<double>, 1> & mesh, const nda::array<std::complex<double>, 1> & data);

    nda::array<double, 1> evaluate(const nda::array<double, 1> & grid);

    private:
    std::vector<complex_t> _phis;
    std::vector<matrix_t> _abcds;

  };

//  typedef solver<mpf_class> SolverMpf;
} // namespace nevanlinna
