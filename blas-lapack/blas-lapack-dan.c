#include <dan.h>
#include <io.h>
#ifdef _HAVE_BLAS_LAPACK
#include <cblas.h>
#endif
#include <assert.h>

/*
  Wrapper functions for BLAS and LAPACK linear algebra routines.
  
  Currently, XTX and matprod malloc the answer themselves; that
  might be poor style, especially if you want the answer to form part
  of a larger matrix. Nevertheless it seems that they need to start with
  zeroed memory, so it makes sure you don't forget that.
  
  This guy's website is very helpful on calling BLAS and LAPACK:
  http://seehuhn.de/pages/linear
*/

#ifdef _HAVE_BLAS_LAPACK
static double dlamch(char CMACH){
    extern double dlamch_(char *CMACHp);
    return dlamch_(&CMACH);
}

void eigen(double *A, int dim, int nvecs, bool values_only, double *evals, double *evecs) {

    /* 
       Compute eigendecomposition of square matrix A.
       A has 'dim' rows and columns and is assumed to be stored in column major order, 
       with the increment in memory to jump between columns also equal to dim (if A were a submatrix of
       a larger matrix, that wouldn't be true and LDA would need to be passed in and set correctly).
       
       Sufficient space for the evecs and evals must be allocated outside this function,
       but all other work arrays are alloc()d and free()d here (including the potentially informative ISUPPZ)
       
       The argument nvecs specifies the number of evecs/evals you want to compute.
       if( nvecs >= dim || nvecs <= 0) then all evecs and evals are computed. 
       Otherwise the first nvecs evecs/evals are computed (considered in order of decreasing eval).

       Note that in the output the evals are in *increasing* order, which might not be what you expect.
       The evecs are in a corresponding order.
    */
    
    extern void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, int *Np,
			double *A, int *LDAp, double *VLp, double *VUp,
			int *ILp, int *IUp, double *ABSTOLp, int *Mp,
			double *W, double *Z, int *LDZp, int *ISUPPZ,
			double *WORK, int *LWORKp, int *IWORK, int *LIWORKp,
			int *INFOp);
    
    char
	JOBZ = values_only ? 'N' : 'V',    // 'N' doesn't compute evecs (faster)
	RANGE,                             // see below
	UPLO = 'L' ;                       // use upper / lower triangle of A ?
    int 
	LDA = dim,                         /* increment between columns (NB if e.g. A is some sort of submatrix
					      of a larger matrix, this might not be the same as dim */
	IL, IU,                            // see below
	M,                                 // will contain total number of evals found
	LDZ = dim,

	*ISUPPZ,                           // indices (2 for each eval) indicating support of evectors
	                                   // (?meaning) 

	LWORK,                             // size of WORK
	*IWORK,                            // work array of ints
	LIWORK,                            // size of IWORK
	INFO,                              // receives error code from fortran routine
	tmp_int ;                          // used to get size of work array of ints
    double
	VL, VU,                            // see below
	ABSTOL=dlamch('S'),                                      
	*WORK,                             // work array of doubles
	tmp_double ;                       // used to get size of work array of ints
    /*
      RANGE
      'A' computes all evals/evecs and ignores IL, IU, VL, VU
      'I' computes ILth through IUth evals/evecs 
          (considering evals ordered in increasing order of eval size; numbering from 1)
      'V' computes evals/evecs whose evals lie within (VL, VU]
    */
    
    if (nvecs > 0 && nvecs < dim) {
	RANGE = 'I' ;
	IL = dim - nvecs + 1;
	IU = dim ;
    }
    else  {
	RANGE = 'A' ;
	nvecs = dim ;
	IL = IU = -1.0 ;
    }
    VL = VU = -1.0 ;
    
    /* First request optimal sizes of work arrays */

    LWORK = LIWORK = -1 ;
    
       dsyevr_(&JOBZ, &RANGE, &UPLO, &dim, A, &LDA, &VL, &VU,
	    &IL, &IU, &ABSTOL, &M, evals, evecs,  &LDZ, ISUPPZ,
	    &tmp_double, &LWORK, &tmp_int, &LIWORK, &INFO);
    if(INFO != 0) {
	PRINT("INFO = %d", INFO) ;
	ERROR("eigen: dsyevr returned with error code when computing optimal size of work arrays.\n") ;
    }

    LWORK = (int) tmp_double ;
    LIWORK = tmp_int ;
    // PRINT("LWORK = %d\nLIWORK = %d\n", LWORK, LIWORK) ;
    
    WORK = ALLOC(LWORK, double) ;
    IWORK = ALLOC(LIWORK, int) ;
    ISUPPZ = ALLOC(2 * nvecs, int) ;
    
    /* Now do real thing */
    
    dsyevr_(&JOBZ, &RANGE, &UPLO, &dim, A, &LDA, &VL, &VU,
	    &IL, &IU, &ABSTOL, &M, evals, evecs,  &LDZ, ISUPPZ,
	    WORK, &LWORK, IWORK, &LIWORK, &INFO);

    if(INFO != 0) {
	PRINT("INFO = %d", INFO) ;
	ERROR("eigen: dsyevr returned with error code.\n") ;
    }
    
    FREE(WORK) ;
    FREE(IWORK) ;
    FREE(ISUPPZ) ;
}

double *XTX(double *X, int nrow, int ncol, double alpha) {
    
    /*
      Compute alpha X'X
      where alpha is a scalar multiplier
    */

    double *ans = CALLOC(ncol * ncol, double) ; if(ans == NULL) ERROR("XTX alloc error\n") ;
    
    cblas_dsyrk(CblasColMajor,
		CblasLower,
		CblasTrans,
		ncol,        // N = 'order' (?) of resulting matrix
		nrow,        // K = nrows
		alpha,         // alpha = scalar multiplier of XTX
		X,           // A = K x N matrix
		nrow,        // lda = nrows
		0.0,         // beta = scalar multiplier of C (don't worry about it)
		ans,
		ncol) ;
    
    return ans ;
}

#endif
double *matprod_general(double *A, double *B, bool Atranspose, bool Btranspose,
			int nrowA, int ncolA, int nrowB, int ncolB, double alpha, double beta, double *C) {

    /*
      Compute alpha LR + beta C
      where alpha, beta are scalar multipliers, and
      L = A' if Atranspose else L = A, and
      R = B' if Btranspose else R = B

      This currently assumes that storage is column-major.
    
      This is more general than matprod and should replace it.
    */

    int nrowL, ncolL,  // dimensions of left and right matrices 
	nrowR, ncolR ; // after any requested transpositions
    int K ;

    if(Atranspose) {
	nrowL = ncolA ;
	ncolL = nrowA ;
    }
    else {
	nrowL = nrowA ;
	ncolL = ncolA ;
    }
    if(Btranspose) {
	nrowR = ncolB ;
	ncolR = nrowB ;
    }
    else {
	nrowR = nrowB ;
	ncolR = ncolB ;
    }

    if(ncolL != nrowR) {
	PRINT("matprod non-conformable error!\nncolL = %d\tnrowR = %d\n", ncolL, nrowR) ;
	exit(1) ;
    }
    K = ncolL ;

#ifdef _HAVE_BLAS_LAPACK
    cblas_dgemm(CblasColMajor,
		Atranspose ? CblasTrans : CblasNoTrans,
		Btranspose ? CblasTrans : CblasNoTrans,
		nrowL, // M
		ncolR, // N
		K,     // The dimension that collapses, i.e. over which dot products are taken
		alpha,   // alpha = scalar multiplier of product
		A,
		nrowA, // LDA = A memory increment from column to column
		B,
		nrowB, // LDB = B memory increment from column to column
		beta,   // beta = multiplier of matrix to add to product
		C,     // C = matrix to add to product (after scaling)
		nrowL  // LDC = memory increment from column to column
		) ;
#else
    /* This is my budget version, for use when linking against BLAS is problematic. */

    int i, j, k, ij ;

    if( !Atranspose && Btranspose ) {
	/* So we're computing AB': ans has dimensions nrowA x nrowB (column major). */
	assert(ncolA = ncolB) ;
    
	for(i = 0 ; i < nrowA ; i++) {
	    for(j = 0 ; j < nrowB ; j++) {
		ij = j * nrowA + i ;
		C[ij] *= beta ;
		for(k = 0 ; k < ncolA ; k++)
		    C[ij] +=  A[k * nrowA + i] * B[k * nrowB + j] ; // A_ik * B_jk ;
		C[ij] *= alpha ;
	    }
	}
    }
    else if( Atranspose && !Btranspose ) {
	/* So we're computing A'B: ans has dimensions ncolA x ncolB (column major). */
	assert(nrowA = nrowB) ;
	
	for(i = 0 ; i < ncolA ; i++) {
	    for(j = 0 ; j < ncolB ; j++) {
		ij = j * ncolA + i ;
		C[ij] *= beta ;
		for(k = 0 ; k < nrowA ; k++)
		    C[ij] +=  A[nrowA * i + k] * B[nrowB * j + k] ; // A_ki * B_kj ;
		C[ij] *= alpha ;
	    }
	}
    }
    else
	ERROR("AB' and A'B have been implemented in Dan's budget matrix product code, but not AB and A'B'. You should not see this message if you are linking against BLAS/LAPACK") ;

    
#endif
    
    return C ;
}


double *matprod(double *A, double *B, bool Atranspose, bool Btranspose,
		int nrowA, int ncolA, int nrowB, int ncolB, double alpha) {

    /*
      Compute alpha LR
      where alpha is a scalar multiplier, and
      L = A' if Atranspose else L = A, and
      R = B' if Btranspose else R = B

      This currently assumes that storage is column-major.
    */

    int nrowL, ncolL,  // dimensions of left and right matrices 
	nrowR, ncolR ; // after any requested transpositions
    int K ;
    double *ans, beta ;

    if(Atranspose) {
	nrowL = ncolA ;
	ncolL = nrowA ;
    }
    else {
	nrowL = nrowA ;
	ncolL = ncolA ;
    }
    if(Btranspose) {
	nrowR = ncolB ;
	ncolR = nrowB ;
    }
    else {
	nrowR = nrowB ;
	ncolR = ncolB ;
    }

    if(ncolL != nrowR) {
	PRINT("matprod non-conformable error!\nncolL = %d\tnrowR = %d\n", ncolL, nrowR) ;
	exit(1) ;
    }
    K = ncolL ;

    ans = CALLOC(nrowL * ncolR, double) ; // C = matrix to add to product (after scaling)
    beta = 0.0 ; // beta = multiplier of matrix to add to product
    
    /* 
       2008-03-XX
       Since I've now created matprod_general (descended from original version of this function) 
       which makes use of more of the functionality of
       cblas_dgemm, I've redefined this as a call to matprod_general. Specifically, matprod_general allows
       a matrix C times scalar beta to be added to the product.
       It also does not allocate the result; the storage for C is used.
    */

    return matprod_general(A, B, Atranspose, Btranspose, nrowA, ncolA, nrowB, ncolB, alpha, beta, ans) ;
    ERROR("matprod now redefined as call to matprod general") ;
    
    /*
    cblas_dgemm(CblasColMajor,
		Atranspose ? CblasTrans : CblasNoTrans,
		Btranspose ? CblasTrans : CblasNoTrans,
		nrowL, // M
		ncolR, // N
		K,     // The dimension that collapses, i.e. over which dot products are taken
		alpha,   // alpha = scalar multiplier of product
		A,
		nrowA, // LDA = A memory increment from column to column
		B,
		nrowB, // LDB = B memory increment from column to column
		beta,   // beta = multiplier of matrix to add to product
		ans,   // C = matrix to add to product (after scaling)
		nrowL  // LDC = memory increment from column to column
		) ;
    */
    return ans ;
}

int main2() {

    double A[] = { // 3 x 2
	1,2,3,
	4,5,6
    } ;

    double B[] = { // 2 x 4
	1,2,
	3,4,
	5,6,
	7,8
    } ;

    double *ans = matprod(A, B, FALSE, FALSE, 3, 2, 2, 4, 1.0) ;
    
    write_matrix_double(ans, stdout, 3, 4, "%4.0lf") ;

    FREE(ans) ;
    return 0 ;
}

#ifdef _HAVE_BLAS_LAPACK
int main1() {
    
    double *ans, *evecs, *evals ;
    
    double X[] = {  // 3 x 4 matrix, column major (e.g 3 SNPs, 4 individs)
	3, 1, 3,
	1, 5, 9,
	2, 6, 5,
	1.3, 2.9, 0.01
    };
    
    double Y[] = {  // 3 x 2 matrix, column major (e.g 3 SNPs, 2 individs)
	0.5, .998, 12.1,
	11.2, 7.3, 8.995
    } ;
    
    double Z[] = {
	1, 2, 3,
	4, 5, 6,
	7, 8, 9
    } ;


    printf("X\n") ;
    write_matrix_double(X, stdout, 3, 4, "%.3lf ") ;
    putchar('\n') ;


    ans = XTX(X, 3, 4, 1.0) ;
    
    printf("XTX\n") ;
    write_matrix_double(ans, stdout, 4, 4, "%.5lf ") ;
    
    FREE(ans) ;

    putchar('\n') ;
    printf("t(X) %%*%% X\n") ;
    ans = matprod(X, X, TRUE, FALSE, 3, 4, 3, 4, 1.0) ;
    write_matrix_double(ans, stdout, 4, 4, "%.5lf ") ;
        
    FREE(ans) ;

    printf("Y\n") ;
    write_matrix_double(Y, stdout, 3, 2, "%.3lf ") ;

    putchar('\n') ;
    printf("t(X) %%*%% Y\n") ;
    ans = matprod(X, Y, TRUE, FALSE, 3, 4, 3, 2, 1.0) ;
    write_matrix_double(ans, stdout, 4, 2, "%.5lf ") ;
    

    printf("Z\n") ;
    write_matrix_double(Z, stdout, 3, 3, "%.3lf ") ;

    evals = ALLOC(3, double) ;
    evecs = ALLOC(3*3, double) ;
    eigen(Z, 3, 3, FALSE, evals, evecs) ;

    PRINT("eigenvalues of Z\n\n") ;
    // write_evals(evals, stdout, 3) ;
    
    PRINT("eigenvectors of Z\n\n") ;
    // write_evecs(evecs, stdout, 3, 3) ;
    
    return 0 ;
}
#endif
