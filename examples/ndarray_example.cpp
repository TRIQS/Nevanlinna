
#include "triqs_Nevanlinna/kernels.hpp"
#include <mpi/mpi.hpp>
#include <nda/traits.hpp>
#include <nda/h5.hpp>

using namespace std::literals;

int main(int argc, char **argv) {
  mpi::environment env(argc, argv);
  triqs_Nevanlinna::Nevanlinna_kernel kernel;

  size_t N_omega = 5000;
  double omega_min = -10;
  double omega_max = 10;
  double eta = 0.1;
  nda::array<std::complex<double>, 3> G_iw;
  nda::array<std::complex<double>, 1> mesh;

  h5::file input(std::string(DATA_PATH) + "/input.h5", 'r');
  h5::group top(input);
  nda::h5_read(top, std::string("data"), G_iw);
  nda::array<double, 1> tmp;
  nda::h5_read(top, std::string("mesh"), tmp);
  mesh = tmp;
  mesh *= 1i;

  kernel.init(mesh, G_iw);

  nda::array<std::complex<double>, 1> grid(N_omega);
  int i = 0;
  std::transform(grid.begin(), grid.end(), grid.begin(), [&](const std::complex<double> & ) {return omega_min + (i++)*(omega_max - omega_min)/(N_omega - 1) + eta*1i;});

  auto G_omega = kernel.evaluate(grid);
}