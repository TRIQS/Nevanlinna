#ifndef TRIQS_NEVANLINNA_TYPES_HPP
#define TRIQS_NEVANLINNA_TYPES_HPP

#include <complex>
#include <Eigen/Dense>
#include <nda/nda.hpp>

#ifdef WITH_MPFR
#include <boost/multiprecision/mpfr.hpp>
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

namespace triqs_Nevanlinna {
    using nda::range;

    static constexpr auto _ = nda::range::all;
    static constexpr auto __ = nda::ellipsis{};

    static constexpr int mp_digits = 100;
#ifdef WITH_MPFR
    using complex_mpt = std::complex<boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<mp_digits>, boost::multiprecision::et_off>>;
#else
    using complex_mpt = std::complex<boost::multiprecision::number<boost::multiprecision::cpp_dec_float<mp_digits>, boost::multiprecision::et_off>>;
#endif
    using matrix_cplx_mpt  = Eigen::Matrix<complex_mpt, Eigen::Dynamic, Eigen::Dynamic>;
}

#endif //TRIQS_NEVANLINNA_TYPES_HPP
