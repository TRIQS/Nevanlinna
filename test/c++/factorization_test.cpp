#include <gtest/gtest.h>
#include <nda/gtest_tools.hpp>
#include <h5/h5.hpp>
#include <nda/h5.hpp>
#include <triqs_Nevanlinna/Nevanlinna_factorization.hpp>

using namespace triqs_Nevanlinna;

TEST(NevanlinnaFactorization, EvaluateData) { // NOLINT
  std::string root = TEST_PATH;
  Nevanlinna_factorization a;
  h5::file data(root + "/data.h5", 'r');
  nda::vector<std::complex<double>> im_data;
  nda::vector<std::complex<double>> im_grid;
  nda::vector<double> re_data;
  nda::vector<double> re_grid;
  h5_read(data, "input/data", im_data);
  h5_read(data, "input/grid", im_grid);
  h5_read(data, "output/data", re_data);
  h5_read(data, "output/grid", re_grid);

  a.build(im_grid, im_data);
  nda::vector<double> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}

TEST(NevanlinnaFactorization, NegativeImGrid) { // NOLINT
  std::string root = TEST_PATH;
  Nevanlinna_factorization a;
  h5::file data(root + "/data.h5", 'r');
  nda::vector<std::complex<double>> im_data;
  nda::vector<std::complex<double>> im_grid;
  h5_read(data, "input/data", im_data);
  h5_read(data, "input/grid", im_grid);
  im_grid *= -1;
  ASSERT_THROW(a.build(im_grid, im_data), Nevanlinna_negative_grid_error);
}

TEST(NevanlinnaFactorization, UnitializedContinuation) { // NOLINT
  std::string root = TEST_PATH;
  Nevanlinna_factorization a;
  h5::file data(root + "/data.h5", 'r');
  nda::vector<std::complex<double>> re_grid;
  h5_read(data, "output/grid", re_grid);
  ASSERT_THROW(a.evaluate(re_grid), Nevanlinna_uninitialized_error);
}

TEST(NevanlinnaFactorization, CheckEta) { // NOLINT
  std::string root = TEST_PATH;
  Nevanlinna_factorization a;
  h5::file data(root + "/data.h5", 'r');
  nda::vector<std::complex<double>> im_data;
  nda::vector<std::complex<double>> im_grid;
  nda::vector<std::complex<double>> re_grid;

  h5_read(data, "input/data", im_data);
  h5_read(data, "input/grid", im_grid);
  h5_read(data, "output/grid", re_grid);
  nda::vector<std::complex<double>> re_data(re_grid.shape(0));
  double eta = 0.5;
  re_grid += eta*1i;
  std::transform(re_grid.begin(), re_grid.end(), re_data.begin(), [&](const std::complex<double> & w) {return 1.0/(w);});

  a.build(im_grid, im_data);
  nda::vector<std::complex<double>> result = a.evaluate(re_grid);
  ASSERT_TRUE(array_are_close(re_data, result, 1e-10));
}
