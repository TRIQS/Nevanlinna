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

  template <typename T>
  concept ArithmeticTypes = std::is_arithmetic<T>::value && std::is_convertible<T, real_mpt>::value;

  template <typename T>
  concept CastableTypes = std::is_convertible<T, real_mpt>::value && !std::is_same<T, real_mpt>::value;

  template <typename T>
  concept ConvertibleTypes = std::is_convertible<T, real_mpt>::value;
} // namespace triqs_Nevanlinna

namespace std {
  // std::complex explicit template specialization for triqs_Nevanlinna::real_mpt to overcome libc++ limitations.
  template <> class complex<triqs_Nevanlinna::real_mpt> {
    using real_t = triqs_Nevanlinna::real_mpt;

    public:
    complex() : _real(real_t(0)), _imag(real_t(0)) {}
    complex(real_t real, real_t imag) : _real(real), _imag(imag) {}
    complex(real_t real) : _real(real), _imag(real_t(0)) {}

    template <triqs_Nevanlinna::ArithmeticTypes S> complex(S real) : _real(real_t(real)), _imag(real_t(0)) {}
    template <triqs_Nevanlinna::ArithmeticTypes S> complex(const std::complex<S> &cplx) : _real(real_t(cplx.real())), _imag(real_t(cplx.imag())) {}

    template <triqs_Nevanlinna::CastableTypes S> complex<real_t> &operator=(const S &rhs) {
      _real = real_t(rhs);
      _imag = real_t(0);
      return *this;
    }
    template <triqs_Nevanlinna::CastableTypes S> complex<real_t> &operator=(const std::complex<S> &rhs) {
      _real = real_t(rhs.real());
      _imag = real_t(rhs.imag());
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator-=(const complex<S> &rhs) {
      _real -= real_t(rhs.real());
      _imag -= real_t(rhs.imag());
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator+=(const complex<S> &rhs) {
      _real += real_t(rhs.real());
      _imag += real_t(rhs.imag());
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator*=(const complex<S> &rhs) {
      *this = *this * rhs;
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator/=(const complex<S> &rhs) {
      *this = *this / rhs;
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator-=(const S &rhs) {
      _real -= real_t(rhs);
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator+=(const S &rhs) {
      _real += real_t(rhs);
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator*=(const S &rhs) {
      _real *= real_t(rhs);
      _imag *= real_t(rhs);
      return *this;
    }

    template <triqs_Nevanlinna::ConvertibleTypes S> complex<real_t> &operator/=(const S &rhs) {
      _real /= real_t(rhs);
      _imag /= real_t(rhs);
      return *this;
    }

    complex<real_t> operator-() const { return complex<real_t>(-_real, -_imag); }

    private:
    real_t _real;
    real_t _imag;

    public:
    real_t real() const { return _real; }
    real_t imag() const { return _imag; }

    void real(const real_t &re) { _real = re; }
    void imag(const real_t &im) { _imag = im; }
  };

  template <> inline bool operator==(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return x.real() == y.real() && x.imag() == y.imag();
  }

  template <triqs_Nevanlinna::CastableTypes S> inline bool operator==(const complex<triqs_Nevanlinna::real_mpt> &x, const S &y) {
    return x.real() == y && x.imag() == 0;
  }

  template <triqs_Nevanlinna::CastableTypes S> inline bool operator==(const S &y, const complex<triqs_Nevanlinna::real_mpt> &x) {
    return x.real() == y && x.imag() == 0;
  }

  inline auto operator*(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    triqs_Nevanlinna::real_mpt a = x.real();
    triqs_Nevanlinna::real_mpt b = x.imag();
    triqs_Nevanlinna::real_mpt c = y.real();
    triqs_Nevanlinna::real_mpt d = y.imag();
    return typename std::complex<triqs_Nevanlinna::real_mpt>((a * c - b * d), (a * d + b * c));
  }

  inline auto operator/(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    triqs_Nevanlinna::real_mpt a     = x.real();
    triqs_Nevanlinna::real_mpt b     = x.imag();
    triqs_Nevanlinna::real_mpt c     = y.real();
    triqs_Nevanlinna::real_mpt d     = y.imag();
    triqs_Nevanlinna::real_mpt denom = c * c + d * d;
    triqs_Nevanlinna::real_mpt r     = (a * c + b * d) / denom;
    triqs_Nevanlinna::real_mpt i     = (b * c - a * d) / denom;
    return typename std::complex<triqs_Nevanlinna::real_mpt>(r, i);
  }

  inline auto operator+(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    triqs_Nevanlinna::real_mpt a = x.real();
    triqs_Nevanlinna::real_mpt b = x.imag();
    triqs_Nevanlinna::real_mpt c = y.real();
    triqs_Nevanlinna::real_mpt d = y.imag();
    return typename std::complex<triqs_Nevanlinna::real_mpt>((a + c), (b + d));
  }

  inline auto operator-(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    triqs_Nevanlinna::real_mpt a = x.real();
    triqs_Nevanlinna::real_mpt b = x.imag();
    triqs_Nevanlinna::real_mpt c = y.real();
    triqs_Nevanlinna::real_mpt d = y.imag();
    return typename std::complex<triqs_Nevanlinna::real_mpt>((a - c), (b - d));
  }

  inline triqs_Nevanlinna::real_mpt abs(const complex<triqs_Nevanlinna::real_mpt> &x) { return sqrt(x.real() * x.real() + x.imag() * x.imag()); }

  inline complex<triqs_Nevanlinna::real_mpt> conj(const complex<triqs_Nevanlinna::real_mpt> &x) {
    return complex<triqs_Nevanlinna::real_mpt>(x.real(), -x.imag());
  }

  inline complex<triqs_Nevanlinna::real_mpt> polar(const triqs_Nevanlinna::real_mpt &rho, const triqs_Nevanlinna::real_mpt &arg) {
    triqs_Nevanlinna::real_mpt x = boost::multiprecision::cos(arg);
    triqs_Nevanlinna::real_mpt y = boost::multiprecision::sin(arg);
    x *= rho;
    y *= rho;
    return complex<triqs_Nevanlinna::real_mpt>(x, y);
  }

  inline triqs_Nevanlinna::real_mpt arg(const complex<triqs_Nevanlinna::real_mpt> &x) {
    triqs_Nevanlinna::real_mpt arg_;
    arg_ = boost::multiprecision::atan2(x.imag(), x.real());
    return arg_;
  }

  inline complex<triqs_Nevanlinna::real_mpt> sqrt(const complex<triqs_Nevanlinna::real_mpt> &x) {
    return std::polar(boost::multiprecision::sqrt(std::abs(x)), std::arg(x) / triqs_Nevanlinna::real_mpt(2));
  }

  inline triqs_Nevanlinna::real_mpt real(const complex<triqs_Nevanlinna::real_mpt> &x) { return x.real(); }

  inline triqs_Nevanlinna::real_mpt imag(const complex<triqs_Nevanlinna::real_mpt> &x) { return x.imag(); }
} // namespace std

#endif //TRIQS_NEVANLINNA_COMPLEX_HPP
