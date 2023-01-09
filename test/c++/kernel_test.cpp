#include <triqs/test_tools/gfs.hpp>
#include <triqs_Nevanlinna/Nevanlinna_kernel.hpp>

using namespace triqs_Nevanlinna;

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
