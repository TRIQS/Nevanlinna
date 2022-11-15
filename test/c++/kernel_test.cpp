#include <triqs/test_tools/gfs.hpp>
#include <triqs_Nevanlinna/kernel.hpp>

using namespace triqs_Nevanlinna;

TEST(NevanlinnaKernel, EvaluateData) { // NOLINT
  std::string root = TEST_PATH;
  kernel a;
  h5::file data(root + "/data.h5", 'r');
  nda::array<std::complex<double>, 1> im_data;
  nda::array<std::complex<double>, 1> im_grid;
  nda::array<double, 1> re_data;
  nda::array<double, 1> re_grid;
  h5_read(data, "input/data", im_data);
  h5_read(data, "input/grid", im_grid);
  h5_read(data, "output/data", re_data);
  h5_read(data, "output/grid", re_grid);

  a.solve(im_grid, im_data);
  nda::array<double, 1> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}

TEST(NevanlinnaKernel, NegativeImGrid) { // NOLINT
  std::string root = TEST_PATH;
  kernel a;
  h5::file data(root + "/data.h5", 'r');
  nda::array<std::complex<double>, 1> im_data;
  nda::array<std::complex<double>, 1> im_grid;
  h5_read(data, "input/data", im_data);
  h5_read(data, "input/grid", im_grid);
  im_grid *= -1;
  ASSERT_THROW(a.solve(im_grid, im_data), Nevanlinna_negative_grid_error);
}

TEST(NevanlinnaKernel, CheckEta) { // NOLINT
  std::string root = TEST_PATH;
  kernel a;
  h5::file data(root + "/data.h5", 'r');
  nda::array<std::complex<double>, 1> im_data;
  nda::array<std::complex<double>, 1> im_grid;
  nda::array<std::complex<double>, 1> re_grid;

  h5_read(data, "input/data", im_data);
  h5_read(data, "input/grid", im_grid);
  h5_read(data, "output/grid", re_grid);
  nda::array<std::complex<double>, 1> re_data(re_grid.shape(0));
  double eta = 0.5;
  re_grid += eta*1i;
  std::transform(re_grid.begin(), re_grid.end(), re_data.begin(), [&](const std::complex<double> & w) {return 1.0/(w);});

  a.solve(im_grid, im_data);
  nda::array<std::complex<double>, 1> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}
