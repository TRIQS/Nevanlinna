(changelog)=

# Changelog

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
