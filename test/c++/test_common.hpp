#ifndef NEVANLINNA_TEST_COMMON
#define NEVANLINNA_TEST_COMMON

#include <nda/mpi.hpp>
#include <nda/gtest_tools.hpp>

#define MAKE_MAIN_MPI                                                                                                                                \
  int main(int argc, char **argv) {                                                                                                                  \
    ::testing::InitGoogleTest(&argc, argv);                                                                                                          \
    if (mpi::has_env) {                                                                                                                              \
      mpi::environment env(argc, argv);                                                                                                              \
      std::cout << "MPI environment detected\n";                                                                                                     \
      return RUN_ALL_TESTS();                                                                                                                        \
    } else                                                                                                                                           \
      return RUN_ALL_TESTS();                                                                                                                        \
  }

#endif // NEVANLINNA_TEST_COMMON
