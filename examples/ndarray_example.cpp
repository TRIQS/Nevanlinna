#include <triqs_Nevanlinna/kernels.hpp>
#include <mpi/mpi.hpp>
#include <h5/h5.hpp>
#include <nda/nda.hpp>
#include <nda/h5.hpp>

using namespace std::literals;

int main(int argc, char **argv) {
  // Initialize MPI environment
  mpi::environment env(argc, argv);

  // Create kernel object
  triqs_Nevanlinna::Nevanlinna_kernel kernel;

  // Define real-frequency grid
  size_t N_w   = 5000;
  double w_min = -10., w_max = 10.;
  double eta = 0.1;
  auto del   = (w_max - w_min) / (N_w - 1);
  auto grid  = nda::basic_array{w_min + nda::arange(N_w) * del + eta * 1i};

  // Define imaginary time grid and data arrays
  auto mesh = nda::array<std::complex<double>, 1>{};
  auto G_iw = nda::array<std::complex<double>, 3>{};

  // Read imaginary time data from "input.h5" file
  auto input = h5::file(DATA_PATH + "/input.h5"s, 'r');
  h5::read(input, "data", G_iw);
  h5::read(input, "mesh", mesh);

  // Build Nevanlinna factorization
  kernel.init(mesh, G_iw);

  // Perform analytical continuation onto real frequency axis
  auto G_w = kernel.evaluate(grid);
}
