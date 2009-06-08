#include <dan.h>
#include <unistd.h>
#include <io.h>
#include <iogeno.h>
#include <assert.h>
#include <blas-lapack-dan.h>


/*
  Use pre-computed SNP loadings to project genotype data into principal component space.
  This code descended from snpload.c, modification by Celine Bellenguez.

  Matrix is A x B means that it has A rows and B columns.  Matrices
  are stored column major, i.e. consecutive values in computer memory
  run down columns.

  n = number of individuals
  L = number of SNP loci
  nvecs = number of eigenvectors

  SNP loadings V is nvecs x L
  genotype data X is n x L

  This code computes XV' which is n x nvecs

  Thus position of individual i on principal component v is element (i, v) of
  XV'
*/


#define VERBOSE 0

int main(int argc, char *argv[]) {
    int c, n_tot, n_inc, nvecs, a, b, L_inc, L_tot, l, i ;
    char *geno_file, *evec_file, *freq_file, *outfile ; //, outfile[250] ;
    double *v, *x, *xx, *y, *freq ;
    bool *snp_include, *indiv_include, *evec_include, rescale=FALSE ;
    FILE *f ;

    while((c = getopt(argc, argv, "a:b:N:v:L:g:e:f:o:s")) != -1) {
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
	    outfile = optarg ; break ;
	case 's':
	    rescale = TRUE ; break ;
	case '?':
	    ERROR("Unrecognised option") ;
	}
    }

    assert(a >= 0 && b >= a) ;

    // PRINT("%s\t/* Construct inclusion vectors */\n", timestring()) ; fflush(stdout) ;
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
        
    PRINT("# %s\tComputing projection of genotype data into Principal Coordinate space\n",
	  timestring()) ; fflush(stdout) ;
    PRINT("# number of SNPs = %d\n", L_tot) ;
    // PRINT("# [a,b) = [%d,%d) -> L_inc = %d\n", a, b, L_inc) ;
    PRINT("# number of individuals = %d\n", n_tot) ;
    // PRINT("n_inc = %d\n", n_inc) ;
    PRINT("# number of PCs = %d\n", nvecs) ;
    PRINT("# genotype data file = %s\n", geno_file) ;
    PRINT("# Principal components file = %s\n", evec_file) ;
    PRINT("# Allele frequency file = %s\n", freq_file) ;
    // PRINT("outfile = %s\n", outfile) ;
    
    if(VERBOSE)
	PRINT("%s\t/* Read eigenvectors */\n", timestring()) ; fflush(stdout) ;
    // v = nvecs x L_inc, column major
    v = ALLOC(L_inc * nvecs, double) ; // was nvecs * ninc in snpload.c
    f = fopen(evec_file, "r") ;
    read_submatrix_double(f, nvecs, evec_include, L_tot, snp_include,  "%lf", v) ; // differs from snpload.c
    fclose(f) ;
    
    if(VERBOSE)
	PRINT("%s\t/* Read genotypes */\n", timestring()) ; fflush(stdout) ;
    // x is n_tot x L_inc, column major
    x = ALLOC(L_inc * n_inc, double) ;
    f = fopen(geno_file, "r") ;
    read_submatrix_genotypes_double(f, n_tot, indiv_include, L_tot, snp_include,  x) ;
    fclose(f) ;

    if(VERBOSE) {
	PRINT("%s\t/* First %d entries of genotypes x: */\n", timestring(), MIN(n_inc*nvecs, 100)) ;
	for(i = 0 ; i < MIN(n_inc*nvecs, 100) ; i++) PRINT("%lf ", x[i]) ;
	putchar('\n') ;
    }

    if(VERBOSE)
	PRINT("%s\t/* Read freqs */\n", timestring()) ; fflush(stdout) ;
    freq = ALLOC(L_inc, double) ;
    f = fopen(freq_file, "r") ;
    read_submatrix_double(f, L_tot, snp_include, 1, NULL, "%lf", freq) ;
    fclose(f) ;
    
    if(VERBOSE)
	PRINT("%s\t/* Set NAs to freq */\n", timestring()) ; fflush(stdout) ;
    xx = x ;
    for(l = 0 ; l < L_inc ; l++)
	for(i = 0 ; i < n_inc ; i++, xx++) {
	    if( rescale ) {
		if(*xx == MISSING) *xx = 0.0 ;
		else {
		    *xx -= 2 * freq[l];
		    *xx /= sqrt(2 * freq[l] * (1 - freq[l]));
		}
	    }
	    else if(*xx == MISSING) *xx = 2*freq[l] ;
	}
    if(VERBOSE) {
	PRINT("%s\t/* First %d entries of freqs: */\n", timestring(), MIN(L_inc, 100)) ;
	for(i = 0 ; i < MIN(L_inc, 100) ; i++) PRINT("%lf ", freq[i]) ;
	putchar('\n') ;
    }

    if(VERBOSE)
	PRINT("%s\t/* Compute projection */\n", timestring()) ; fflush(stdout) ;

    /* Now compute XV' which is n_tot x nvecs */
    y = matprod(x, v, FALSE, TRUE, n_tot, L_inc, nvecs, L_inc, 1.0) ; // differs from snpload.c

    // sprintf(outfile, "%s-%d-%d", outfile_partial, a+1, b) ;
    f = fopen(outfile, "w") ;
    write_matrix_double(y, f, n_tot, nvecs, "%lf ") ; // differs from snpload.c
    fclose(f) ;
    
    FREE(y) ;
    FREE(snp_include) ;
    FREE(v) ;
    FREE(x) ;
    FREE(freq) ;
    PRINT("%s\t/* Done */\n", timestring()) ;

    return 0 ;
}
