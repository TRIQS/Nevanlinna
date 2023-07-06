#include <triqs/test_tools/gfs.hpp>
#include <triqs_Nevanlinna/types.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <Eigen/Dense>

using namespace triqs_Nevanlinna;

TEST(NevanlinnaComplexMath, SimpleMath) {
  std::complex<real_mpt> x;
  std::complex<real_mpt> y;
  x = std::complex(2.0, 1.0);
  y = 2;
  auto z = x+y;
  ASSERT_NEAR(z.real().convert_to<double>(), 4, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 1, 1e-12);
  z = x-y;
  ASSERT_NEAR(z.real().convert_to<double>(), 0, 1e-12);
  z = x*y;
  ASSERT_NEAR(z.real().convert_to<double>(), 4, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 2, 1e-12);
  z = x/y;
  ASSERT_NEAR(z.real().convert_to<double>(), 1, 1e-12);
}

TEST(NevanlinnaComplexMath, MixedMath) {
  std::complex<real_mpt> x;
  std::complex<double> y;
  double w;
  x = std::complex(2.0, 1.0);
  y = std::complex(2.0, 2.0);
  w = 2.0;
  std::complex<real_mpt> z = x + y;
  ASSERT_NEAR(z.real().convert_to<double>(), 4, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 3, 1e-12);
  z = x-y;
  ASSERT_NEAR(z.real().convert_to<double>(), 0, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), -1, 1e-12);
  z = x*y;
  ASSERT_NEAR(z.real().convert_to<double>(), 2, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 6, 1e-12);
  z = x/y;
  ASSERT_NEAR(z.real().convert_to<double>(), 0.75, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), -0.25, 1e-12);

  z = y+x;
  ASSERT_NEAR(z.real().convert_to<double>(), 4, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 3, 1e-12);
  z = y-x;
  ASSERT_NEAR(z.real().convert_to<double>(), 0, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 1, 1e-12);
  z = y*x;
  ASSERT_NEAR(z.real().convert_to<double>(), 2, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 6, 1e-12);
  z = y/x;
  ASSERT_NEAR(z.real().convert_to<double>(), 1.2, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 0.4, 1e-12);

  z = w+x;
  ASSERT_NEAR(z.real().convert_to<double>(), 4, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 1, 1e-12);
  z = w-x;
  ASSERT_NEAR(z.real().convert_to<double>(), 0, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), -1, 1e-12);
  z = w*x;
  ASSERT_NEAR(z.real().convert_to<double>(), 4, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), 2, 1e-12);
  z = w/x;
  ASSERT_NEAR(z.real().convert_to<double>(), 0.8, 1e-12);
  ASSERT_NEAR(z.imag().convert_to<double>(), -0.4, 1e-12);
}

TEST(NevanlinnaComplexMath, InplaceMath) {
  std::complex<real_mpt> x;
  std::complex<real_mpt> y;
  std::complex<real_mpt> z;
  x = std::complex(2.0, 1.0);
  y = std::complex(2.0, 2.0);
  x/=y;
  ASSERT_NEAR(x.real().convert_to<double>(), 0.75, 1e-12);
  ASSERT_NEAR(x.imag().convert_to<double>(), -0.25, 1e-12);
  x = std::complex(2.0, 1.0);
  z = std::complex(2.0, 2.0);
  x/=z;
  ASSERT_NEAR(x.real().convert_to<double>(), 0.75, 1e-12);
  ASSERT_NEAR(x.imag().convert_to<double>(), -0.25, 1e-12);
}

TEST(NevanlinnaComplexMath, MatrixMath) {
  using complex_t = std::complex<real_mpt>;
  using matrix_t = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
  using namespace triqs_Nevanlinna;

  matrix_t x(2,2);
  matrix_t y(2,2);
  x << 1., 2., 2., 1.;
  y << -1.0/3.0, 2.0/3.0, 2.0/3.0, -1.0/3.0;
  matrix_t z = x.inverse();
  ASSERT_NEAR(z(0,0).real().convert_to<double>(), y(0,0).real().convert_to<double>(), 1e-9);
}

TEST(NevanlinnaComplexMath, EigenSolverMath) {
  using complex_t = std::complex<real_mpt>;
  using matrix_t = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
  using namespace triqs_Nevanlinna;
  using namespace std::complex_literals;

  matrix_t M(2,2);
  M << 1., 2., 1.i, 1.;
  Eigen::ComplexEigenSolver<matrix_t> ces;
  ces.compute(M);
  matrix_t D = ces.eigenvalues();
  ASSERT_NEAR(D(1,0).real().convert_to<double>(), 2, 1e-9);
  ASSERT_NEAR(D(1,0).imag().convert_to<double>(), 1, 1e-9);
}
