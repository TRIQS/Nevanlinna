#include <triqs_Nevanlinna/kernels.hpp>
#include <mpi/mpi.hpp>
#include <h5/h5.hpp>
#include <nda/nda.hpp>
#include <nda/h5.hpp>

using namespace std::complex_literals;

int main(int argc, char **argv) {
  // Initialize MPI environment
  mpi::environment env(argc, argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Create kernel object
  triqs_Nevanlinna::Caratheodory_kernel kernel;

  // Define real-frequency grid
  size_t N_w   = 500;
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
  mesh *= 1.i;
  // Build Caratheodory factorization
  kernel.init(mesh, G_iw);
  double end = MPI_Wtime();
  if(!rank) std::cout<<"Nevan init is done in "<<end -start<<std::endl;
  // Perform analytical continuation onto real frequency axis
  auto G_w = kernel.evaluate(grid);
  start = MPI_Wtime();
  if(!rank) std::cout<<"Nevan eval is done in "<<start - end<<std::endl;

  // Create kernel object
  triqs_Nevanlinna::Caratheodory_kernel cakernel;
  N_w   = 500;
  del   = (w_max - w_min) / (N_w - 1);
  grid  = nda::basic_array{w_min + nda::arange(N_w) * del + eta * 1i};
  start = MPI_Wtime();
  cakernel.init(mesh, G_iw);
  end = MPI_Wtime();
  if(!rank) std::cout<<"Cara init is done in "<<end -start<<std::endl;
  // Perform analytical continuation onto real frequency axis
  auto G_ij_w = cakernel.evaluate(grid);
  start = MPI_Wtime();
  if(!rank) std::cout<<"Cara eval is done in "<<start - end<<std::endl;
}
