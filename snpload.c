#include <dan.h>
#include <unistd.h>
#include <io.h>
#include <iogeno.h>
#include <assert.h>
#include <blas-lapack-dan.h>

/*
  Read eigenvectors and compute loadings for SNPs in specified range [a,b).

  Matrix is A x B means that it has A rows and B columns.  Matrices are
  stored column major, i.e. consecutive values in computer memory run
  down columns.

  n = number of individuals
  L = number of SNP loci
  nvecs = number of eigenvectors

  genotype data X is n x L
  eigenvectors V is n x nvecs 

  This code computes V'X which is nvecs x L

  Thus loading of SNP l on principal component v is element (v, l) of
  V'X.
*/

int main(int argc, char *argv[]) {
    int c, n_tot, n_inc, nvecs, a, b, L_inc, L_tot, l, i ;
    char *geno_file, *evec_file, *freq_file, *outfile_partial, outfile[250] ;
    double *v, *x, *xx, *y, *freq ;
    bool *snp_include, *indiv_include, *evec_include ;
    FILE *f ;

    while((c = getopt(argc, argv, "a:b:N:v:L:g:e:f:o:")) != -1) {
	switch(c) {
	case 'a':
	    a = atoi(optarg) - 1; break ;
	case 'b':
	    b = atoi(optarg) ; break ;
	case 'N':
	    n_tot = atoi(optarg) ; break ;
	case 'v':
	    nvecs = atoi(optarg) ; break ;
	case 'L':
	    L_tot = atoi(optarg) ; break ;
	case 'g':
	    geno_file = optarg ; break ;
	case 'e':
	    evec_file = optarg ; break ;
	case 'f':
	    freq_file = optarg ; break ;
	case 'o':
	    outfile_partial = optarg ; break ;
	case '?':
	    ERROR("Unrecognised option") ;
	}
    }

    assert(a >= 0 && b >= a) ;

    PRINT("%s\t/* Construct inclusion vectors */\n", timestring()) ; fflush(stdout) ;
    indiv_include = NULL ;
    n_inc = n_tot ;
    evec_include = NULL ;
    snp_include = ALLOC(L_tot, bool) ;
    L_inc = 0 ;
    for(l = 0 ; l < L_tot ; l++) {
	snp_include[l] = (l >= a && l < b) ;
	if( snp_include[l] ) L_inc++ ;
    }
    assert(L_inc == b - a) ;
        
    PRINT("L_tot = %d\n", L_tot) ;
    PRINT("[a,b) = [%d,%d) -> L_inc = %d\n", a, b, L_inc) ;
    PRINT("n_tot = %d\n", n_tot) ;
    PRINT("n_inc = %d\n", n_inc) ;
    PRINT("nvecs = %d\n", nvecs) ;
    PRINT("geno file = %s\n", geno_file) ;
    PRINT("evec file = %s\n", evec_file) ;
    PRINT("freq file = %s\n", freq_file) ;
    PRINT("outfile_partial = %s\n", outfile_partial) ;
    
    PRINT("%s\t/* Read eigenvectors */\n", timestring()) ; fflush(stdout) ;
    // v is n_tot x nvecs, column major
    v = ALLOC(nvecs * n_inc, double) ;
    f = fopen(evec_file, "r") ;

    read_submatrix_double(f, n_tot, indiv_include, nvecs, evec_include, "%lf", v) ;

    fclose(f) ;
    
    // PRINT("%s\t/* First %d entries of eigenvectors v: */\n", timestring(), MIN(n_inc*nvecs, 100)) ;
    // for(i = 0 ; i < MIN(n_inc*nvecs, 100) ; i++) PRINT("%lf ", v[i]) ;
    // putchar('\n') ;
    
    PRINT("%s\t/* Read genotypes */\n", timestring()) ; fflush(stdout) ;
    // x is n_tot x L_inc, column major
    x = ALLOC(L_inc * n_inc, double) ;
    f = fopen(geno_file, "r") ;
    read_submatrix_genotypes_double(f, n_tot, indiv_include, L_tot, snp_include,  x) ;
    fclose(f) ;

    // PRINT("%s\t/* First %d entries of genotypes x: */\n", timestring(), MIN(n_inc*nvecs, 100)) ;
    // for(i = 0 ; i < MIN(n_inc*nvecs, 100) ; i++) PRINT("%lf ", x[i]) ;
    // putchar('\n') ;
    
    PRINT("%s\t/* Read freqs */\n", timestring()) ; fflush(stdout) ;
    freq = ALLOC(L_inc, double) ;
    f = fopen(freq_file, "r") ;
    read_submatrix_double(f, L_tot, snp_include, 1, NULL, "%lf", freq) ;
    fclose(f) ;
    
    PRINT("%s\t/* Set NAs to freq */\n", timestring()) ; fflush(stdout) ;
    xx = x ;
    for(l = 0 ; l < L_inc ; l++)
	for(i = 0 ; i < n_inc ; i++)
	    if(*xx++ == MISSING) *(xx - 1) = 2 * freq[l] ;
    
    PRINT("%s\t/* First %d entries of freqs: */\n", timestring(), MIN(L_inc, 100)) ;
    for(i = 0 ; i < MIN(L_inc, 100) ; i++) PRINT("%lf ", freq[i]) ;
    putchar('\n') ;

    PRINT("%s\t/* Compute SNP loadings */\n", timestring()) ; fflush(stdout) ;
    y = matprod(v, x, TRUE, FALSE, n_tot, nvecs, n_tot, L_inc, 1.0) ;
    sprintf(outfile, "%s-%d-%d", outfile_partial, a+1, b) ;
    f = fopen(outfile, "w") ;
    write_matrix_double(y, f, nvecs, L_inc, "%lf\t") ;
    fclose(f) ;
    
    FREE(y) ;
    FREE(snp_include) ;
    FREE(v) ;
    FREE(x) ;
    FREE(freq) ;
    PRINT("%s\t/* Done! */\n", timestring()) ;

    return 0 ;
}
