#ifndef TRIQS_NEVANLINNA_COMPLEX_HPP
#define TRIQS_NEVANLINNA_COMPLEX_HPP

#ifdef WITH_MPFR
#include <boost/multiprecision/mpfr.hpp>
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

namespace triqs_Nevanlinna {
  static constexpr int mp_digits = 100;
#ifdef WITH_MPFR
  using real_mpt = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_off>;
#else
  using real_mpt = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<mp_digits>, boost::multiprecision::et_off>;
#endif
} // namespace triqs_Nevanlinna

namespace std {

  // std::complex explicit template specialization for triqs_Nevanlinna::real_mpt to overcome libc++ limitations.
  template <> class complex<triqs_Nevanlinna::real_mpt> {

    public:
    using value_type = triqs_Nevanlinna::real_mpt;

    complex() = default;
    complex(value_type re) : _re(re) {}
    complex(value_type re, value_type im) : _re(re), _im(im) {}

    template <typename S>
      requires(is_arithmetic_v<S>)
    complex(S re) : _re(re) {}

    template <typename S>
      requires(is_arithmetic_v<S>)
    complex(const complex<S> &cplx) : _re(cplx.real()), _im(cplx.imag()) {}

    template <convertible_to<value_type> S> complex &operator=(const complex<S> &rhs) {
      _re = rhs.real();
      _im = rhs.imag();
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator-=(const complex<S> &rhs) {
      _re -= rhs.real();
      _im -= rhs.imag();
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator+=(const complex<S> &rhs) {
      _re += rhs.real();
      _im += rhs.imag();
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator*=(const complex<S> &rhs) {
      *this = *this * rhs;
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator/=(const complex<S> &rhs) {
      *this = *this / rhs;
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator-=(const S &rhs) {
      _re -= rhs;
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator+=(const S &rhs) {
      _re += rhs;
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator*=(const S &rhs) {
      _re *= rhs;
      _im *= rhs;
      return *this;
    }

    template <convertible_to<value_type> S> complex &operator/=(const S &rhs) {
      _re /= rhs;
      _im /= rhs;
      return *this;
    }

    complex operator-() const { return complex(-_re, -_im); }

    friend inline auto operator*(const complex &x, const complex &y) {
      auto [xr, xi] = x;
      auto [yr, yi] = y;
      return complex((xr * yr - xi * yi), (xr * yi + xi * yr));
    }

    friend inline auto operator/(const complex &x, const complex &y) {
      auto [xr, xi] = x;
      auto [yr, yi] = y;
      auto denom    = yr * yr + yi * yi;
      auto r        = (xr * yr + xi * yi) / denom;
      auto i        = (xi * yr - xr * yi) / denom;
      return complex(r, i);
    }

    friend inline auto operator+(const complex &x, const complex &y) {
      auto [xr, xi] = x;
      auto [yr, yi] = y;
      return complex((xr + yr), (xi + yi));
    }

    friend inline auto operator-(const complex &x, const complex &y) {
      auto [xr, xi] = x;
      auto [yr, yi] = y;
      return complex((xr - yr), (xi - yi));
    }

    friend inline bool operator==(const complex &x, const complex &y) = default;

    template <convertible_to<value_type> S> friend inline bool operator==(const complex &x, const S &y) { return x.real() == y && x.imag() == 0; }

    template <convertible_to<value_type> S> friend inline bool operator==(const S &y, const complex &x) { return x.real() == y && x.imag() == 0; }

    value_type real() const { return _re; }
    value_type imag() const { return _im; }

    friend inline auto real(const complex &x) { return x._re; }
    friend inline auto imag(const complex &x) { return x._im; }

    private:
    value_type _re = 0;
    value_type _im = 0;
  };

  inline auto abs(const complex<triqs_Nevanlinna::real_mpt> &x) { return sqrt(x.real() * x.real() + x.imag() * x.imag()); }
  inline auto conj(const complex<triqs_Nevanlinna::real_mpt> &x) { return complex<triqs_Nevanlinna::real_mpt>(x.real(), -x.imag()); }

  inline auto polar(const triqs_Nevanlinna::real_mpt &rho, const triqs_Nevanlinna::real_mpt &arg) {
    triqs_Nevanlinna::real_mpt x = boost::multiprecision::cos(arg);
    triqs_Nevanlinna::real_mpt y = boost::multiprecision::sin(arg);
    x *= rho;
    y *= rho;
    return complex<triqs_Nevanlinna::real_mpt>(x, y);
  }

  inline auto arg(const complex<triqs_Nevanlinna::real_mpt> &x) {
    triqs_Nevanlinna::real_mpt arg_;
    arg_ = boost::multiprecision::atan2(x.imag(), x.real());
    return arg_;
  }

  inline auto sqrt(const complex<triqs_Nevanlinna::real_mpt> &x) {
    return polar(boost::multiprecision::sqrt(abs(x)), arg(x) / triqs_Nevanlinna::real_mpt(2));
  }

} // namespace std

#endif //TRIQS_NEVANLINNA_COMPLEX_HPP
