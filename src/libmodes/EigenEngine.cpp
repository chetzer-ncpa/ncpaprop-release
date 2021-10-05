#include "slepceps.h"
#include "slepcst.h"
#include <complex>
#include <string>
#include <iostream>
#include "EigenEngine.h"



int NCPA::EigenEngine::doESSCalculation( double *diag, int Nz_grid, double dz, double tol, 
	double *k_min, double *k_max, PetscInt *nconv, double *k2, double **v ) {

	//
	// Declarations related to Slepc computations
	//
	Mat            A;           // problem matrix
	EPS            eps;         // eigenproblem solver context
	ST             stx;
	KSP            kspx;
	PC             pcx;
	EPSType        type;        // CHH 191022: removed const qualifier
	PetscReal      re, im;
	PetscScalar    kr, ki, *xr_;
	Vec            xr, xi;
	PetscInt       Istart, Iend, col[3], its, maxit;
	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar    value[3];	
	PetscErrorCode ierr;
	PetscMPIInt    rank, size;
	int i, j, nev;
	double h2 = dz * dz;


	// Initialize Slepc
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);  

	// Create the matrix A to use in the eigensystem problem: Ak=kx
	ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,Nz_grid,Nz_grid); CHKERRQ(ierr);
	ierr = MatSetFromOptions(A); CHKERRQ(ierr);

	// the following Preallocation call is needed in PETSc version 3.3
	ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL); CHKERRQ(ierr);
	// or use: ierr = MatSetUp(A); 

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Compute the operator matrix that defines the eigensystem, Ax=kx
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	// Make matrix A 
	ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
	if (Istart==0) 
		FirstBlock=PETSC_TRUE;
	if (Iend==Nz_grid) 		/* @todo check if should be Nz_grid-1 */
		LastBlock=PETSC_TRUE;    

	value[0]=1.0/h2; 
	value[2]=1.0/h2;
	for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++ ) {
		value[1] = -2.0/h2 + diag[i];
		col[0]=i-1; 
		col[1]=i; 
		col[2]=i+1;
		ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES); CHKERRQ(ierr);
	}
	if (LastBlock) {
		i=Nz_grid-1; 
		col[0]=Nz_grid-2; 
		col[1]=Nz_grid-1;
		ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
	}
	if (FirstBlock) {
		i=0; 
		col[0]=0; 
		col[1]=1; 
		value[0]=-2.0/h2 + diag[0]; 
		value[1]=1.0/h2;
		ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
	}

	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// CHH 191022: MatGetVecs() is deprecated, changed to MatCreateVecs()
	ierr = MatCreateVecs(A,PETSC_NULL,&xr); CHKERRQ(ierr);
	ierr = MatCreateVecs(A,PETSC_NULL,&xi); CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Create the eigensolver and set various options
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* 
	Create eigensolver context
	*/
	ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);

	/* 
	Set operators. In this case, it is a standard eigenvalue problem
	*/
	ierr = EPSSetOperators(eps,A,PETSC_NULL); CHKERRQ(ierr);
	ierr = EPSSetProblemType(eps,EPS_HEP); CHKERRQ(ierr);

	/*
	Set solver parameters at runtime
	*/
	ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
	ierr = EPSSetType(eps,"krylovschur"); CHKERRQ(ierr);
	ierr = EPSSetDimensions(eps,10,PETSC_DECIDE,PETSC_DECIDE); CHKERRQ(ierr); // leaving this line in speeds up the code; better if this is computed in chunks of 10? - consult Slepc manual
	ierr = EPSSetTolerances(eps,tol,PETSC_DECIDE); CHKERRQ(ierr);

	ierr = EPSGetST(eps,&stx); CHKERRQ(ierr);
	ierr = STGetKSP(stx,&kspx); CHKERRQ(ierr);
	ierr = KSPGetPC(kspx,&pcx); CHKERRQ(ierr);
	ierr = STSetType(stx,"sinvert"); CHKERRQ(ierr);
	ierr = KSPSetType(kspx,"preonly");
	ierr = PCSetType(pcx,"cholesky");
	ierr = EPSSetInterval(eps,pow(*k_min,2),pow(*k_max,2)); CHKERRQ(ierr);
	ierr = EPSSetWhichEigenpairs(eps,EPS_ALL); CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Solve the eigensystem
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = EPSSolve(eps);CHKERRQ(ierr);
	/*
	Optional: Get some information from the solver and display it
	*/
	ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
	ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
	ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
	ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	Display solution and clean up
	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* 
	Get number of converged approximate eigenpairs
	*/
	ierr = EPSGetConverged(eps,nconv);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %d\n\n",*nconv);CHKERRQ(ierr);

	
	if ((*nconv)>0) {
		for (i=0;i<(*nconv);i++) {
			ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
#if defined(PETSC_USE_COMPLEX)
			re = PetscRealPart(kr);
			im = PetscImaginaryPart(kr);
#else
			re = kr;
			im = ki;
#endif 
			k2[(*nconv)-i-1] = re; // proper count of modes
			ierr = VecGetArray(xr,&xr_);CHKERRQ(ierr);
			for (j = 0; j < Nz_grid; j++) {
				v[j][(*nconv)-i-1] = xr_[j]/sqrt(dz); //per Slepc the 2-norm of xr_ is=1; we need sum(v^2)*dz=1 hence the scaling xr_/sqrt(dz)
			}
			ierr = VecRestoreArray(xr,&xr_);CHKERRQ(ierr);
		}
	}

	ierr = EPSDestroy(&eps);CHKERRQ(ierr);
	ierr = MatDestroy(&A);  CHKERRQ(ierr);
	ierr = VecDestroy(&xr); CHKERRQ(ierr);
	ierr = VecDestroy(&xi); CHKERRQ(ierr); 

	

	return 0;
}

int NCPA::EigenEngine::doWideAngleCalculation( int Nz_grid, double dz, double k_min, double k_max,
			double tol, int nev, double *kd, double *md, double *cd, 
			PetscInt *nconv, double *kH, double **v, std::string disp_msg ) {

	Mat            A, B;       
	EPS            eps;  		// eigenproblem solver context      
	ST             stx;
	EPSType  type;		// CHH 191028: removed const qualifier
	PetscReal      re, im;
	PetscScalar    kr, ki, *xr_;
	Vec            xr, xi;
	PetscInt       Istart, Iend, col[3], its, maxit;
	PetscBool      FirstBlock=PETSC_FALSE, LastBlock=PETSC_FALSE;
	PetscScalar    value[3];
	PetscErrorCode ierr;
	PetscMPIInt    rank, size;
	Vec            V_SEQ;
	VecScatter     ctx;

	int    i, j;
	double h2 = dz * dz;
	double sigma = 0.5*(k_min+k_max);
    int nev_2 = nev*2;
    int n_2   = Nz_grid*2;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

    if (rank == 0) {
        std::cout << disp_msg << std::endl;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    //   Compute the operator matrices that define the eigensystem, A*x = k.B*x
    //   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n_2,n_2);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
    ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,n_2,n_2);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);CHKERRQ(ierr);
    
    // the following Preallocation calls are needed in PETSc version 3.3
    ierr = MatSeqAIJSetPreallocation(A, 3, PETSC_NULL);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(B, 2, PETSC_NULL);CHKERRQ(ierr);

    /*
    We solve the quadratic eigenvalue problem (Mk^2 + Ck +D)v = 0 where k are the 
    eigenvalues and v are the eigenvectors and M , C, A are N by N matrices.
    We linearize this by denoting u = kv and arrive at the following generalized
    eigenvalue problem:

     / -D  0 \  (v )      / C  M \  (v ) 
    |         | (  ) = k |        | (  )
     \ 0   M /  (kv)      \ M  0 /  (kv)
     
     ie. of the form
     A*x = k.B*x
     so now we double the dimensions of the matrices involved: 2N by 2N.
     
     Matrix D is tridiagonal and has the form
     Main diagonal     : -2/h^2 + omega^2/c(1:N)^2 + F(1:N)
     Upper diagonal    :  1/h^2
     Lower diagonal    :  1/h^2
     Boundary condition:  A(1,1) =  (1/(1+h*beta) - 2)/h^2 + omega^2/c(1)^2 + F(1)
     
     where 
     F = 1/2*rho_0"/rho_0 - 3/4*(rho_0')^2/rho_0^2 where rho_0 is the ambient
     stratified air density; the prime and double prime means first and second 
     derivative with respect to z.
     beta  = alpha - 1/2*rho_0'/rho_0
     alpha is given in
     Psi' = alpha*Psi |at z=0 (i.e. at the ground). Psi is the normal mode.
     If the ground is rigid then alpha is zero. 
     
     Matrix M is diagonal:
     Main diagonal  : u0(1:N)^2/c(1:N)^2 - 1
      u0 = scalar product of the wind velocity and the horizontal wavenumber k_H
      u0 = v0.k_H
      
     Matrix C is diagonal:
     Main diagonal: 2*omega*u0(1:N)/c(1:N)^2 
     
    */

    // Assemble the A matrix (2N by 2N)
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
    if (Istart==0) FirstBlock=PETSC_TRUE;
    if (Iend==n_2) LastBlock =PETSC_TRUE;

    // matrix -D is placed in the first N by N block
    // kd[i]=(omega/c_T)^2
    for( i=(FirstBlock? Istart+1: Istart); i<(LastBlock? (Iend/2)-1: Iend/2); i++ ) {
        value[0]=-1.0/h2; value[1] = 2.0/h2 - kd[i]; value[2]=-1.0/h2;
        col[0]=i-1; col[1]=i; col[2]=i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (LastBlock) {
        i=(n_2/2)-1; col[0]=(n_2/2)-2; col[1]=(n_2/2)-1; value[0]=-1.0/h2; value[1]=2.0/h2 - kd[(n_2/2)-1];
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    // boundary condition
    if (FirstBlock) {
        i=0; col[0]=0; col[1]=1; value[0]=2.0/h2 - kd[0]; value[1]=-1.0/h2;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    // Insert matrix M into the lower N by N block of A
    // md is u0^2/c^2-1
    for ( i=(n_2/2); i<n_2; i++ ) {
        ierr = MatSetValue(A,i,i,md[i-(n_2/2)],INSERT_VALUES); CHKERRQ(ierr); 
    }

    // Assemble the B matrix
    for ( i=0; i<(n_2/2); i++ ) {
        col[0]=i; col[1]=i+(n_2/2); value[0]=cd[i]; value[1]=md[i];
        ierr = MatSetValues(B,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    for ( i=(n_2/2); i<n_2; i++ ) {
        ierr = MatSetValue(B,i,i-(n_2/2),md[i-(n_2/2)],INSERT_VALUES); CHKERRQ(ierr); 
    }

    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    // CHH 191028: MatGetVecs is deprecated, changed to MatCreateVecs
    ierr = MatCreateVecs(A,PETSC_NULL,&xr);CHKERRQ(ierr);
    ierr = MatCreateVecs(A,PETSC_NULL,&xi);CHKERRQ(ierr);
    //ierr = MatGetVecs(A,PETSC_NULL,&xr);CHKERRQ(ierr);
    //ierr = MatGetVecs(A,PETSC_NULL,&xi);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  Create the eigensolver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Create eigensolver context
    */
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

    /* 
       Set operators. In this case, it is a quadratic eigenvalue problem
    */
    ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps,EPS_GNHEP);CHKERRQ(ierr);

    /*
       Set solver parameters at runtime
    */
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
    ierr = EPSSetType(eps,"krylovschur"); CHKERRQ(ierr);
    ierr = EPSSetDimensions(eps,nev_2,PETSC_DECIDE,PETSC_DECIDE); CHKERRQ(ierr);
    ierr = EPSSetTarget(eps,sigma); CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps,tol,PETSC_DECIDE); CHKERRQ(ierr);

    ierr = EPSGetST(eps,&stx); CHKERRQ(ierr);
    ierr = STSetType(stx,"sinvert"); CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE); CHKERRQ(ierr);
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        Solve the eigensystem
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = EPSSolve(eps);CHKERRQ(ierr);
    /*
       Optional: Get some information from the solver and display it
    */
    ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
    ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps,&nev_2,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev_2);CHKERRQ(ierr);
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Display solution and clean up
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /* 
       Get number of converged approximate eigenpairs
    */
    ierr = EPSGetConverged(eps,nconv); CHKERRQ(ierr);
    
    if ((*nconv)>0) {
        for( i=0; i<(*nconv); i++ ) {
            ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
            #ifdef PETSC_USE_COMPLEX
                re = PetscRealPart(kr);
                im = PetscImaginaryPart(kr);
            #else
                re = kr;
                im = ki;
            #endif 
            kH[i] = re;
            ierr = VecScatterCreateToAll(xr,&ctx,&V_SEQ);CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx,xr,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx,xr,V_SEQ,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
            if (rank == 0) {
                ierr = VecGetArray(V_SEQ,&xr_);CHKERRQ(ierr);
                for (j = 0; j < Nz_grid; j++) {
                    v[j][i] = xr_[j]/sqrt(dz);
                }
                ierr = VecRestoreArray(V_SEQ,&xr_);CHKERRQ(ierr);
            }
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            ierr = VecDestroy(&V_SEQ);CHKERRQ(ierr);
        }
    }

    // Free work space
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A);  CHKERRQ(ierr);
    ierr = MatDestroy(&B);  CHKERRQ(ierr);
    ierr = VecDestroy(&xr); CHKERRQ(ierr);
    ierr = VecDestroy(&xi); CHKERRQ(ierr); 

    return 0;

}