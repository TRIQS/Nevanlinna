file(GLOB_RECURSE sources *.cpp)

if(NOT TRIQS_Support)
    list(FILTER sources EXCLUDE REGEX ".*solver_core\\.cpp")
endif()

add_library(${PROJECT_NAME}_c ${sources})
add_library(${PROJECT_NAME}::${PROJECT_NAME}_c ALIAS ${PROJECT_NAME}_c)

# Link against nda / triqs and enable warnings
target_link_libraries(${PROJECT_NAME}_c PUBLIC nda::nda_c $<$<BOOL:${TRIQS_Support}>:triqs> PRIVATE $<BUILD_INTERFACE:${PROJECT_NAME}_warnings>)

# Configure target and compilation
set_target_properties(${PROJECT_NAME}_c PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)
target_include_directories(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_include_directories(${PROJECT_NAME}_c SYSTEM INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
target_compile_definitions(${PROJECT_NAME}_c PUBLIC
				NEVANLINNA_GIT_HASH=${PROJECT_GIT_HASH}
				TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:NEVANLINNA_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			  )

# Install library and headers
install(TARGETS ${PROJECT_NAME}_c EXPORT ${PROJECT_NAME}-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")

# ========= Library Dependencies ==========

find_package(Eigen3 REQUIRED CONFIG)
target_link_libraries(${PROJECT_NAME}_c PUBLIC Eigen3::Eigen)
# Prior to Eigen 3.4, the template Eigen::internal::result_of attempts to use
# std::result_of which was removed in C++20 in favor of
# std::invoke_result. However, Eigen also has its own implementation of
# result_of which can be used by simply defining the macro
# EIGEN_HAS_STD_RESULT_OF=0. It doesn't interfere with newer versions of
# Eigen. These will still use std::invoke_result because guarded by yet another
# feature macro.
# https://gitlab.com/libeigen/eigen/-/commit/a31effc3bca2f0924752caeebfd6f61f7edf9a43
target_compile_definitions(${PROJECT_NAME}_c PUBLIC EIGEN_HAS_STD_RESULT_OF=0)

find_package(Boost REQUIRED)
target_link_libraries(${PROJECT_NAME}_c PUBLIC Boost::boost)
# There is a usage of std::unary_function (remove in C++17) in Boost
# container_hash that is not correctly guarded before Boost 1.73.0, because
# _HAS_AUTO_PTR_ETC only exists in MSVC (as far as I can tell).
# https://github.com/boostorg/container_hash/commit/90a0e3663875973909f77062dc274081bf32dc0d
target_compile_definitions(${PROJECT_NAME}_c PUBLIC _HAS_AUTO_PTR_ETC=0)

find_package(OpenMP REQUIRED COMPONENTS CXX)
target_link_libraries(${PROJECT_NAME}_c PUBLIC OpenMP::OpenMP_CXX)

option(WITH_MPFR "Use MPFR as a backend for multiprecision operations." ON)
if(WITH_MPFR)
  find_package(MPFR REQUIRED)
  target_compile_definitions(${PROJECT_NAME}_c PUBLIC WITH_MPFR)
  target_link_libraries(${PROJECT_NAME}_c PUBLIC MPFR)
  install(TARGETS MPFR EXPORT ${PROJECT_NAME}-targets)
endif()

# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(${PROJECT_NAME}_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET ${PROJECT_NAME}_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++23
	  --template=gcc
	  --verbose
	  --force
	  --quiet
          ${sources}
	WORKING_DIRECTORY
	  ${CMAKE_CURRENT_SOURCE_DIR}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()

# ========= Dynamic Analyzer Checks ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
if(ASAN)
  if(NOT TARGET asan)
    find_package(sanitizer REQUIRED COMPONENTS asan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:asan>)
endif()

option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")
if(UBSAN)
  if(NOT TARGET ubsan)
    find_package(sanitizer REQUIRED COMPONENTS ubsan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:ubsan>)
endif()

option(MSAN OFF "Compile library and executables with LLVM Memory Sanitizer")
if(MSAN)
  if(NOT TARGET msan)
    find_package(sanitizer REQUIRED COMPONENTS msan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:msan>)
endif()

option(TSAN OFF "Compile library and executables with LLVM Thread Sanitizer")
if(TSAN)
  if(NOT TARGET tsan)
    find_package(sanitizer REQUIRED COMPONENTS tsan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:tsan>)
endif()
