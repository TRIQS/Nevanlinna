.. highlight:: bash

.. _install:

Install Nevanlinna
******************


Packaged Versions of Nevanlinna
==========================
.. _ubuntu_debian:
Ubuntu Debian packages
----------------------

We provide a Debian package for the Ubuntu LTS Version 22.04 (jammy), which can be installed by following the steps outlined :ref:`here <triqslibs:ubuntu_debian>`, and the subsequent command::

        sudo apt-get install -y triqs_nevanlinna


Compiling Nevanlinna from source
================================

.. note:: To guarantee reproducibility in scientific calculations we strongly recommend the use of a stable `release <https://github.com/TRIQS/triqs/releases>`_ of both TRIQS and its applications.

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/Nevanlinna`` repository from GitHub::

     $ git clone https://github.com/TRIQS/Nevanlinna Nevanlinna.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir Nevanlinna.build && cd Nevanlinna.build

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../Nevanlinna.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Version compatibility
---------------------

Keep in mind that the version of ``Nevanlinna`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd Nevanlinna.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``Nevanlinna`` can be configured using CMake-options::

    cmake ../Nevanlinna.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_Nevanlinna     |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the without TRIQS Support                                 | -DTRIQS_Support=OFF                           |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the without Python Support                                | -DPythonSupport=OFF                           |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the without MPFR                                          | -DWITH_MPFR=OFF                               |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
