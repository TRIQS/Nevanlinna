#ifndef TRIQS_NEVANLINNA_TYPES_HPP
#define TRIQS_NEVANLINNA_TYPES_HPP

#include <Eigen/Dense>
#include <nda/nda.hpp>

#include "triqs_Nevanlinna/complex.hpp"

namespace triqs_Nevanlinna {
  using nda::range;

  static constexpr auto _  = nda::range::all;
  static constexpr auto __ = nda::ellipsis{};

  static constexpr int mp_digits = 100;
#ifdef WITH_MPFR
  using real_mpt = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_off>;
#else
  using real_mpt = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<mp_digits>, boost::multiprecision::et_off>;
#endif
  using complex_mpt     = std::complex<real_mpt>;
  using matrix_cplx_mpt = Eigen::Matrix<complex_mpt, Eigen::Dynamic, Eigen::Dynamic>;
  static const auto I   = complex_mpt{0., 1.};
  static const auto One = complex_mpt{1., 0.};
} // namespace triqs_Nevanlinna
#endif //TRIQS_NEVANLINNA_TYPES_HPP
