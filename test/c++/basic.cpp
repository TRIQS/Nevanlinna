#include <triqs/test_tools/gfs.hpp>
#include <nevanlinna/nevanlinna.hpp>

using namespace nevanlinna;

TEST(Nevanlinna, EvaluateEmpty) { // NOLINT

  nevanlinna_parameters_t p;
  solver a(p);

  ASSERT_ANY_THROW(a.evaluate(nda::array<std::complex<double>, 1>(10)));
}

TEST(Nevanlinna, EvaluateData) { // NOLINT
  std::string root = TEST_PATH;
  nevanlinna_parameters_t p;
  solver a(p);
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


