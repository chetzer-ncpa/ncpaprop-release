# ncpaprop

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5562712.svg)](https://doi.org/10.5281/zenodo.5562712)

# Installation

Run ./configure with appropriate parameters.  Most parameters involve how to link to the [PETSc](https://petsc.org/release/) suite and its [SLEPc](https://slepc.upv.es/) extension.  If, as is likely, you do not already have these both installed with the correct configuration flags (see below), you should use:

	./configure --with-localpetsc --enable-autodependencies
	
Note that the --enable-autodependencies flag is only supported on Linux systems that use apt or yum as package managers.  This will download, configure, and build PETSc and SLEPc locally, within the ncpaprop directory tree (in the extern directory), as part of the configuration process.  Two architectures will be built:

	arch-${OS}-c-real:     --with-scalar-type=real
	arch-${OS}-c-complex:  --with-scalar-type=complex

If you already have architectures of both PETSc and SLEPc using each of these two architectures, you can link to them instead of building them locally.  To do this, set the PETSC_DIR and SLEPC_DIR variables to the root directories of PETSc and SLEPc, and set the PETSC_ARCH_REAL and PETSC_ARCH_COMPLEX variables to the names of the architectures as they were built.  Then, you may use:

	./configure --enable-autodependencies
	
or, alternately, you may specify the values of those four variables in the configure command as:

	./configure PETSC_DIR=${PETSC_DIR} SLEPC_DIR=${SLEP_DIR} PETSC_ARCH_REAL=${PETSC_ARCH_REAL} PETSC_ARCH_COMPLEX=${PETSC_ARCH_COMPLEX} --enable-autodependencies

See the [manual](./NCPA_prop_manual.pdf) for detailed information on additional parameters.

Once the configuration process is complete, simply run 

	make
