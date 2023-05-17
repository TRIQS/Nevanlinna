#include <gtest/gtest.h>
#include <nda/gtest_tools.hpp>
#include <triqs_Nevanlinna/Nevanlinna_kernel.hpp>
#include <triqs_Nevanlinna/Caratheodory_kernel.hpp>

using namespace triqs_Nevanlinna;
using namespace std::complex_literals;

TEST(NevanlinnaKernelFactorization, DiagonalOnly) {
  const double eta = 0.1;
  const double mu1 = 1.0;
  const double mu2 = 1.0;
  const double beta = 10.0;
  const double omega_min = -5.0;
  const double omega_max = 5.0;
  const int n_iw = 100;
  const int n_omega = 100;
  Nevanlinna_kernel a;
  nda::array<std::complex<double>, 2> im_data(n_iw, 2);
  nda::vector<std::complex<double>> im_grid(n_iw);
  nda::array<std::complex<double>, 2> re_data(n_omega, 2);
  nda::vector<std::complex<double>> re_grid(n_omega);

  for(int iw=0, w = -n_iw/2; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2*iw + 1) * M_PI * 1i/beta;
    im_data(iw, 0) = 1./(im_grid(iw) - mu1);
    im_data(iw, 1) = 1./(im_grid(iw) - mu2);
  }

  for(size_t iw = 0; iw < n_omega; ++iw) {
    re_grid(iw) = std::complex<double>(iw*(omega_max - omega_min)/n_omega, eta);
    re_data(iw, 0) = 1./(re_grid(iw) - mu1);
    re_data(iw, 1) = 1./(re_grid(iw) - mu2);
  }

  a.init_diagonal(im_grid, im_data);
  nda::array<std::complex<double>, 2> result = a.evaluate_diagonal(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}
TEST(NevanlinnaKernelFactorization, EvaluateData) {
  const double eta = 0.1;
  const double mu1 = 1.0;
  const double mu2 = 1.0;
  const double beta = 10.0;
  const double omega_min = -5.0;
  const double omega_max = 5.0;
  const int n_iw = 100;
  const int n_omega = 100;
  Nevanlinna_kernel a;
  nda::array<std::complex<double>, 3> im_data(n_iw, 2, 2);
  nda::vector<std::complex<double>> im_grid(n_iw);
  nda::array<std::complex<double>, 3> re_data(n_omega, 2, 2);
  nda::vector<std::complex<double>> re_grid(n_omega);

  for(int iw=0, w = -n_iw/2; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2*iw + 1) * M_PI * 1i/beta;
    im_data(iw, 0, 0) = 1./(im_grid(iw) - mu1);
    im_data(iw, 1, 1) = 1./(im_grid(iw) - mu2);
  }

  for(size_t iw = 0; iw < n_omega; ++iw) {
    re_grid(iw) = std::complex<double>(iw*(omega_max - omega_min)/n_omega, eta);
    re_data(iw, 0, 0) = 1./(re_grid(iw) - mu1);
    re_data(iw, 1, 1) = 1./(re_grid(iw) - mu2);
  }

  a.init(im_grid, im_data);
  nda::array<std::complex<double>, 3> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}

TEST(CaratheodoryKernelFactorization, EvaluateData) {
  const double eta = 0.1;
  const double mu1 = 1.0;
  const double mu2 = -1.0;
  const double beta = 10.0;
  const double omega_min = -5.0;
  const double omega_max = 5.0;
  const int n_iw = 30;
  const int n_omega = 100;
  Caratheodory_kernel a;
  nda::array<std::complex<double>, 3> im_data(n_iw, 2, 2);
  nda::array<std::complex<double>, 1> im_grid(n_iw);
  nda::array<std::complex<double>, 3> re_data(n_omega, 2, 2);
  nda::array<std::complex<double>, 1> re_grid(n_omega);

  for(int iw=0, w = 0; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2*iw + 1) * M_PI * 1i/beta;
    im_data(iw, 0, 0) = 1./(im_grid(iw) - mu1);
    im_data(iw, 1, 1) = 1./(im_grid(iw) - mu2);
  }

  for(size_t iw = 0; iw < n_omega; ++iw) {
    re_grid(iw) = std::complex<double>(omega_min + iw*(omega_max - omega_min)/(n_omega-1), eta);
    re_data(iw, 0, 0) = 1./(re_grid(iw) - mu1);
    re_data(iw, 1, 1) = 1./(re_grid(iw) - mu2);
  }
  a.init(im_grid, im_data);
  nda::array<std::complex<double>, 3> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}

TEST(CaratheodoryKernelFactorization, EvaluateOffDiagonalData) {
  const double eta = 0.1;
  const double beta = 10.0;
  const double omega_min = -5.0;
  const double omega_max = 5.0;
  const int n_iw = 30;
  const int n_omega = 100;
  Eigen::MatrixXcd id = Eigen::MatrixXcd::Identity(2,2);
  Eigen::MatrixXcd H0(2,2);
  H0 << 1.0, 0.5, 0.5, -1.0;
  Caratheodory_kernel a;
  nda::array<std::complex<double>, 3> im_data(n_iw, 2, 2);
  nda::array<std::complex<double>, 1> im_grid(n_iw);
  nda::array<std::complex<double>, 3> re_data(n_omega, 2, 2);
  nda::array<std::complex<double>, 1> re_grid(n_omega);

  for(int iw=0, w = 0; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2*iw + 1) * M_PI * 1i/beta;
    Eigen::MatrixXcd val = (id * im_grid(iw) - H0).inverse();
    im_data(iw, 0, 0) = val(0,0);
    im_data(iw, 0, 1) = val(0,1);
    im_data(iw, 1, 0) = val(1,0);
    im_data(iw, 1, 1) = val(1,1);
  }

  for(size_t iw = 0; iw < n_omega; ++iw) {
    re_grid(iw) = std::complex<double>(omega_min + iw*(omega_max - omega_min)/(n_omega-1), eta);
    Eigen::MatrixXcd val = (id * re_grid(iw) - H0).inverse();
    re_data(iw, 0, 0) = val(0,0);
    re_data(iw, 0, 1) = val(0,1);
    re_data(iw, 1, 0) = val(1,0);
    re_data(iw, 1, 1) = val(1,1);
  }
  a.init(im_grid, im_data);
  nda::array<std::complex<double>, 3> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}

TEST(CaratheodoryKernelFactorization, NonNormalized) {
  const double eta = 0.1;
  const double beta = 10.0;
  const double omega_min = -5.0;
  const double omega_max = 5.0;
  const int n_iw = 30;
  const int n_omega = 100;
  const double norm = 5.0;
  Eigen::MatrixXcd id = Eigen::MatrixXcd::Identity(2,2);
  Eigen::MatrixXcd H0(2,2);
  H0 << 1.0, 0.5, 0.5, -1.0;
  Caratheodory_kernel a;
  nda::array<std::complex<double>, 3> im_data(n_iw, 2, 2);
  nda::array<std::complex<double>, 1> im_grid(n_iw);
  nda::array<std::complex<double>, 3> re_data(n_omega, 2, 2);
  nda::array<std::complex<double>, 1> re_grid(n_omega);

  for(int iw=0, w = 0; iw < n_iw; ++iw, ++w) {
    im_grid(iw) = (2*iw + 1) * M_PI * 1i/beta;
    Eigen::MatrixXcd val = norm*(id * im_grid(iw) - H0).inverse();
    im_data(iw, 0, 0) = val(0,0);
    im_data(iw, 0, 1) = val(0,1);
    im_data(iw, 1, 0) = val(1,0);
    im_data(iw, 1, 1) = val(1,1);
  }

  for(size_t iw = 0; iw < n_omega; ++iw) {
    re_grid(iw) = std::complex<double>(omega_min + iw*(omega_max - omega_min)/(n_omega-1), eta);
    Eigen::MatrixXcd val = norm*(id * re_grid(iw) - H0).inverse();
    re_data(iw, 0, 0) = val(0,0);
    re_data(iw, 0, 1) = val(0,1);
    re_data(iw, 1, 0) = val(1,0);
    re_data(iw, 1, 1) = val(1,1);
  }
  a.init(im_grid, im_data);
  nda::array<std::complex<double>, 3> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}
