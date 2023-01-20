#ifndef NCPAPROP_PETSC_H_INCLUDED
#define NCPAPROP_PETSC_H_INCLUDED

#include "petscvec.h"
#include "petscmat.h"

namespace NCPA {
	void outputVec( Vec &v, double *z, size_t n,
		const std::string &filename ) const;
	void outputSparseMat( Mat &m, size_t nrows,
		const std::string &filename ) const;
}


#endif
