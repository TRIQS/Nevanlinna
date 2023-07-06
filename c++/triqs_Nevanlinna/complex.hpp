#ifndef TRIQS_NEVANLINNA_COMPLEX_HPP
#define TRIQS_NEVANLINNA_COMPLEX_HPP

#ifdef WITH_MPFR
#include <boost/multiprecision/mpfr.hpp>
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

namespace triqs_Nevanlinna {
#ifdef WITH_MPFR
  using real_mpt = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>, boost::multiprecision::et_off>;
#else
  using real_mpt = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<mp_digits>, boost::multiprecision::et_off>;
#endif

  template <typename T>
  concept ArithmeticTypes = std::is_arithmetic<T>::value && std::is_convertible<T, real_mpt>::value;

  template <typename T>
  concept CastableTypes = std::is_convertible<T, real_mpt>::value && !std::is_same<T, real_mpt>::value && std::is_floating_point<T>::value;

  template <typename T>
  concept ConvertibleTypes = std::is_convertible<T, real_mpt>::value;

  template <typename T, typename U>
  concept CompatibleTypes = std::is_convertible<U, T>::value || std::is_convertible<T, U>::value;

  namespace internal {
    template <typename T, CompatibleTypes<T> U> struct complex_return_t {
      using type         = decltype(T{} + U{});
      using complex_type = std::complex<type>;
    };

    template <typename U, typename V>
    inline typename complex_return_t<U, V>::complex_type complex_multiply(const std::complex<U> &lhs, const std::complex<V> &rhs) {
      typename complex_return_t<U, V>::type a = lhs.real();
      typename complex_return_t<U, V>::type b = lhs.imag();
      typename complex_return_t<U, V>::type c = rhs.real();
      typename complex_return_t<U, V>::type d = rhs.imag();
      return typename complex_return_t<U, V>::complex_type((a * c - b * d), (a * d + b * c));
    }

    template <typename U, typename V>
    inline typename complex_return_t<U, V>::complex_type complex_divide(const std::complex<U> &lhs, const std::complex<V> &rhs) {
      typename complex_return_t<U, V>::type a     = lhs.real();
      typename complex_return_t<U, V>::type b     = lhs.imag();
      typename complex_return_t<U, V>::type c     = rhs.real();
      typename complex_return_t<U, V>::type d     = rhs.imag();
      typename complex_return_t<U, V>::type denom = c * c + d * d;
      typename complex_return_t<U, V>::type x     = (a * c + b * d) / denom;
      typename complex_return_t<U, V>::type y     = (b * c - a * d) / denom;
      return typename complex_return_t<U, V>::complex_type(x, y);
    }

    template <typename U, typename V>
    inline typename complex_return_t<U, V>::complex_type complex_add(const std::complex<U> &lhs, const std::complex<V> &rhs) {
      typename complex_return_t<U, V>::type a = lhs.real();
      typename complex_return_t<U, V>::type b = lhs.imag();
      typename complex_return_t<U, V>::type c = rhs.real();
      typename complex_return_t<U, V>::type d = rhs.imag();
      return typename complex_return_t<U, V>::complex_type((a + c), (b + d));
    }

    template <typename U, typename V>
    inline typename complex_return_t<U, V>::complex_type complex_sub(const std::complex<U> &lhs, const std::complex<V> &rhs) {
      typename complex_return_t<U, V>::type a = lhs.real();
      typename complex_return_t<U, V>::type b = lhs.imag();
      typename complex_return_t<U, V>::type c = rhs.real();
      typename complex_return_t<U, V>::type d = rhs.imag();
      return typename complex_return_t<U, V>::complex_type((a - c), (b - d));
    }

    template <typename U, typename V> inline typename complex_return_t<U, V>::complex_type complex_multiply(const std::complex<U> &lhs, V rhs) {
      typename complex_return_t<U, V>::type a = lhs.real();
      typename complex_return_t<U, V>::type b = lhs.imag();
      typename complex_return_t<U, V>::type c = rhs;
      return typename complex_return_t<U, V>::complex_type(a * c, b * c);
    }

    template <typename U, typename V> inline typename complex_return_t<U, V>::complex_type complex_divide(const std::complex<U> &lhs, V rhs) {
      typename complex_return_t<U, V>::type a = lhs.real();
      typename complex_return_t<U, V>::type b = lhs.imag();
      typename complex_return_t<U, V>::type c = rhs;
      typename complex_return_t<U, V>::type x = a / c;
      typename complex_return_t<U, V>::type y = b / c;
      return typename complex_return_t<U, V>::complex_type(x, y);
    }

    template <typename U, typename V> inline typename complex_return_t<U, V>::complex_type complex_divide(V lhs, const std::complex<U> &rhs) {
      typename complex_return_t<U, V>::type a     = lhs;
      typename complex_return_t<U, V>::type b     = 0;
      typename complex_return_t<U, V>::type c     = rhs.real();
      typename complex_return_t<U, V>::type d     = rhs.imag();
      typename complex_return_t<U, V>::type denom = c * c + d * d;
      typename complex_return_t<U, V>::type x     = (a * c) / denom;
      typename complex_return_t<U, V>::type y     = -(a * d) / denom;
      return typename complex_return_t<U, V>::complex_type(x, y);
    }

    template <typename U, typename V> inline typename complex_return_t<U, V>::complex_type complex_add(const std::complex<U> &lhs, V rhs) {
      typename complex_return_t<U, V>::type a = lhs.real();
      typename complex_return_t<U, V>::type b = lhs.imag();
      typename complex_return_t<U, V>::type c = rhs;
      return typename complex_return_t<U, V>::complex_type(a + c, b);
    }

    template <typename U, typename V> inline typename complex_return_t<U, V>::complex_type complex_sub(const std::complex<U> &lhs, V rhs) {
      typename complex_return_t<U, V>::type a = lhs.real();
      typename complex_return_t<U, V>::type b = lhs.imag();
      typename complex_return_t<U, V>::type c = rhs;
      return typename complex_return_t<U, V>::complex_type((a - c), b);
    }
  } // namespace internal
} // namespace triqs_Nevanlinna

namespace std {
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

  template <triqs_Nevanlinna::CastableTypes S> inline auto operator*(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<S> &y) {
    return triqs_Nevanlinna::internal::complex_multiply(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator/(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<S> &y) {
    return triqs_Nevanlinna::internal::complex_divide(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator+(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<S> &y) {
    return triqs_Nevanlinna::internal::complex_add(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator-(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<S> &y) {
    return triqs_Nevanlinna::internal::complex_sub(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator*(const complex<triqs_Nevanlinna::real_mpt> &x, S y) {
    return triqs_Nevanlinna::internal::complex_multiply(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator/(const complex<triqs_Nevanlinna::real_mpt> &x, S y) {
    return triqs_Nevanlinna::internal::complex_divide(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator+(const complex<triqs_Nevanlinna::real_mpt> &x, S y) {
    return triqs_Nevanlinna::internal::complex_add(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator-(const complex<triqs_Nevanlinna::real_mpt> &x, S y) {
    return triqs_Nevanlinna::internal::complex_sub(x, y);
  }
  inline auto operator*(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_multiply(x, y);
  }
  inline auto operator/(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_divide(x, y);
  }
  inline auto operator+(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_add(x, y);
  }
  inline auto operator-(const complex<triqs_Nevanlinna::real_mpt> &x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_sub(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator*(S x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_multiply(y, x);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator/(S x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_divide(x, y);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator+(S x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_add(y, x);
  }
  template <triqs_Nevanlinna::CastableTypes S> inline auto operator-(S x, const complex<triqs_Nevanlinna::real_mpt> &y) {
    return triqs_Nevanlinna::internal::complex_add(-y, x);
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
