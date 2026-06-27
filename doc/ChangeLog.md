# Changelog

## Version 4.0.0

Nevanlinna version 4.0.0 is a compatibility release for TRIQS version 4.0.0.
It is based on the latest app4triqs skeleton and ports the Python bindings from
cpp2py to clair + c2py.

We thank all contributors: Jennifer Coulter, Thomas Hahn, Alexander Hampel, Henri Menke, Dylan Simon, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* Port to TRIQS 4.0 (run the `port_to_triqs4` script; the `triqs.gf` module is now `triqs.gfs`)
* Migrate the Python bindings from cpp2py to clair + c2py, generating the `solver_core` and `kernels_core` modules via `c2py_add_module`
* Use the latest app4triqs skeleton
* Fix the Hardy optimization function
* Fix a compilation error with LLVM 22 on macOS

### cmake
* Use the `clair_c2py_generate_bindings` helper from c2py and disable C++20 module scanning for clair-c2py compatibility
* Suppress `-Wc2y-extensions` warnings (`__COUNTER__` in `SCOPE_EXIT`) for clang

### jenkins
* Migrate to the new k8s-based Jenkins system (#26)
* Install libeigen3-dev and libmpfr-dev in jenkins/Dockerfile

### doc
* Update the doxygen documentation in the kernel and solver headers and bump the Doxyfile to v1.17.0
* Update the tutorials for TRIQS 3.3.x and refresh the landing page
* Add the Flatiron Institute support notice to README.md

## Version 3.3.0

Nevanlinna version 3.3.0 is a compatibility release for TRIQS version 3.3.0

We thank all contributors: Thomas Hahn, Alexander Hampel, Sergei Iskakov, Henri Menke, Nils Wentzell

Find below an itemized list of changes in this release.

### General
* Update install.rst instruction for deb package installation
* fix: add compatibility for older Eigen and Boost
* fix: require Eigen 3.4

### jenkins
* Disable intel build due to boost version issue
* Specify MPFR_ROOT for jenkins osx builds

### cmake
* Undo renaming of cmake config files being installed
* Rename cmake config files being installed
* Raise boost requirement to 1.73

## Version 3.2.0

This is the initial release of Nevanlinna,
an analytic continuation package.
It is compatible against TRIQS Version 3.2.

We thank all contributors: Alexander Hampel, Dylan Simon, Nils Wentzell, Sergei Iskakov
