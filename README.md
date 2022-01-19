# ncpaprop

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5562712.svg)](https://doi.org/10.5281/zenodo.5562712)

**ncpaprop** is a software package aiming at providing a comprehensive set of tested and validated numerical models for simulating the long range propagation of infrasonic signals through the earth’s atmosphere. The algorithms implemented in ncpaprop are designed for frequencies large enough that the effects of buoyancy can be neglected and small enough that propagation to ranges of hundreds to thousands of kilometers is possible without significant signal attenuation. Nominally, **ncpaprop** can, without modification, be used to efficiently model narrowband propagation from 0.1 to 10 Hz and broadband propagation from 0.05 Hz to 2 or 3 Hz. The models become increasingly inefficient with increasing frequency so that run times can become prohibitive if higher frequencies are considered.

The intent behind **ncpaprop** is to provide reliable software engines for the modeling of infrasound propagation rather than a user-friendly working environment. As such no graphical interfaces are included. Rather **ncpaprop** provides a suite of UNIX style command line programs that can be scripted into usersupplied graphical interfaces as desired.

# Prerequisites

The following are required to install **ncpaprop**:

* Bash and Perl interpreters.
* C and C++ compilers.  **ncpaprop** was written and built using the GCC compiler suite.
* GNU Make or an equivalent that supplies the ``make`` command.
* [FFTW](https://www.fftw.org/) for Fourier transforms.
* The [GSL (GNU Scientific Library)](https://www.gnu.org/software/gsl/) for interpolation.

Most Linux distributions should come with the first three; MacOS users may have to install these using XCode and/or an external package manager such as Homebrew.  If using a package manager for **FFTW** and **GSL**, note that the development versions are required.  If you are building **ncpaprop** on a Linux system, you can add the ``--enable-autodependencies`` flag to use the system package manager to install the correct versions of these; the flag supports both **apt** and **yum**.  This flag is not supported on MacOS due to the number of potential package managers to choose from.

# Installation

First, download the repository into the current directory with:

	git clone https://github.com/chetzer-ncpa/ncpaprop-release.git .

Run ``./configure`` with appropriate parameters.  Most parameters involve how to link to the [PETSc](https://petsc.org/release/) suite and its [SLEPc](https://slepc.upv.es/) extension.  If, as is likely, you do not already have these both installed with the correct configuration flags (see below), you should use:

	./configure --with-localpetsc
	
This will download, configure, and build **PETSc** and **SLEPc** locally, within the **ncpaprop** directory tree (in the ``extern`` directory), as part of the configuration process.  Two architectures will be built:

	arch-${OS}-c-real:     --with-scalar-type=real
	arch-${OS}-c-complex:  --with-scalar-type=complex

If you already have architectures of both **PETSc** and **SLEPc** using each of these two architectures, you can link to them instead of building them locally.  To do this, set the ``PETSC_DIR`` and ``SLEPC_DIR`` variables to the root directories of **PETSc** and **SLEPc**, and set the ``PETSC_ARCH_REAL`` and ``PETSC_ARCH_COMPLEX`` variables to the names of the architectures as they were built.  Then, you may use:

	./configure
	
or, alternately, you may specify the values of those four variables in the configure command as:

	./configure PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEP_DIR} PETSC_ARCH_REAL=${PETSC_ARCH_REAL} PETSC_ARCH_COMPLEX=${PETSC_ARCH_COMPLEX}

See the [manual](./NCPA_prop_manual.pdf) for detailed information on additional parameters.

Once the configuration process is complete, simply run 

	make
