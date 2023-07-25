#include <triqs/test_tools/gfs.hpp>
#include <triqs_Nevanlinna/utils.hpp>

using namespace triqs_Nevanlinna;

TEST(test_utils, get_env) {
   int val = get_env_int("NEVANLINNA_TEST_ENV", 1);
   ASSERT_EQ(val, 1);
   setenv("NEVANLINNA_TEST_ENV", "10", 1);
   val = get_env_int("NEVANLINNA_TEST_ENV", 1);
   ASSERT_EQ(val, 10);
};