#include "ncpaprop_petsc.h"
#include <iostream>
#include <fstream>
#include <complex>

#include "petscvec.h"
#include "petscmat.h"

void NCPA::outputVec( Vec &v, double *z, size_t n, const std::string &filename ) {
	PetscScalar *array;
	std::ofstream out( filename );
	out.precision( 12 );
	VecGetArray(v,&array);
	for (size_t i = 0; i < n; i++) {
		if (z != NULL) {
			out << z[i] << "  ";
		}
#ifdef PETSC_USE_COMPLEX
		out << array[i].real() << "  " << array[i].imag() << std::endl;
#else
		out << array[i] << std::endl;
#endif
	}
	VecRestoreArray(v,&array);
	out.close();
}

void NCPA::outputSparseMat( Mat &m, size_t nrows, const std::string &filename ) {
	PetscInt ncols;
	const PetscInt *cols;
	const PetscScalar *vals;
	std::ofstream out( filename );
	out << nrows << std::endl;
	for (size_t i = 0; i < nrows; i++) {
		MatGetRow( m, (PetscInt)i, &ncols, &cols, &vals );
		for (PetscInt j = 0; j < ncols; j++) {
			out 	<< i << " " << cols[j] << " "
#ifdef PETSC_USE_COMPLEX
					<< vals[j].real() << " " << vals[j].imag()
#else
					<< vals[j]
#endif
				<< std::endl;
		}
		MatRestoreRow( m, (PetscInt)i, &ncols, &cols, &vals );
	}
	out.close();
}
