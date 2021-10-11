# ncpaprop

[![DOI](https://zenodo.org/badge/413896280.svg)](https://zenodo.org/badge/latestdoi/413896280)

# Installation

1. Run ./configure with appropriate parameters.  Examples:

To link to an existing PETSc/SLEPc installation:

	./configure PETSC_DIR=/code/petsc SLEPC_DIR=/code/slepc PETSC_ARCH_REAL=arch-linux-c-real PETSC_ARCH_COMPLEX=arch-linux-c-complex --enable-autodependencies

To download and install PETSc and SLEPc locally to the ncpaprop installation:

	./configure --with-localpetsc --enable-autodependencies

See the [manual](./NCPA_prop_manual.pdf) for detailed information on additional parameters.

2. Run 

	make
