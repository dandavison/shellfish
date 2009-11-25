double *matprod_general(double *A, double *B, bool Atranspose, bool Btranspose,
			int nrowA, int ncolA, int nrowB, int ncolB, double alpha, double beta, double *C) ;
double *matprod(double *A, double *B, bool Atranspose, bool Btranspose,
		int nrowA, int ncolA, int nrowB, int ncolB, double alpha) ;
void eigen(double *A, int dim, int nvecs, bool values_only, double *evals, double *evecs) ;
double *XTX(double *X, int nrow, int ncol, double alpha) ;
