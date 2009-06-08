#include <dan.h>
#include <unistd.h>
#include <assert.h>
#include <io.h>
#include <blas-lapack-dan.h>
#include <iogeno.h>

int main(int argc, char *argv[]) {
 
    /*
      Read in data for two sets of individuals, corresponding to one contiguous block of
      the full individual-individual covariance matrix, and compute that block of the
      covariance matrix, after the Patterson et al. centering, zeroing missing data and
      scaling procedure.
    
      Write each column of the sub-matrix (only those elements in lower triangle of full matrix) 
      to a separate file, with a suffix indicating the position of that column in the full column
      of the full cov matrix, so that cat * outputs the full lower triangle of the cov matrix in
      column major order when all submatrices have been computed.
    */
    
    int a, b, c, d, i, i_tot, j_tot, nx, ny, n_tot, n_inc, L_tot, L_inc, ch, nmono ;
    char *genotypes_file=NULL, *data_dir=NULL, *indiv_exclude_file=NULL, *snp_exclude_file=NULL, 
	*freqs_file=NULL, *outdir=NULL, 
	**fnames, outfile[100], *version_string="0.0" ;
    bool symmetric, *indiv_include, *snp_include, *xinclude, *yinclude, verbose=FALSE ;
    FILE *fp ;
    double *x, *y, *freqs_inc, *cov ;

    while((ch = getopt(argc, argv, "a:b:c:d:N:n:L:l:g:i:f:I:S:o:v")) != -1) {
	switch(ch) {
	case 'a':
	    a = atoi(optarg) - 1 ; break ; // indexing from 1 in shell
	case 'b':
	    b = atoi(optarg)     ; break ; // <= in shell
	case 'c':
	    c = atoi(optarg) - 1 ; break ; // indexing from 1 in shell
	case 'd':
	    d = atoi(optarg)     ; break ; // <= 1 in shell
	    //	case 'x':
	    // xchunk = atoi(optarg) ; break ; // not used in C; just for filename
	case 'N':
	    n_tot = atoi(optarg) ; break ;
	    //case 'n':
	    // n_inc = atoi(optarg) ; break ;
	case 'L':
	    L_tot = atoi(optarg) ; break ;
	case 'l':
	    L_inc = atoi(optarg) ; break ;
	case 'g':
	    genotypes_file = optarg ; break ;
	case 'i':
	    data_dir = optarg ; break ;
	case 'f':
	    freqs_file = optarg ; break ;
	case 'I':
	    indiv_exclude_file = optarg ; break ;
	case 'S':
	    snp_exclude_file = optarg ; break ;
	case 'o':
	    outdir = optarg ; break ;
	case 'v':
	    verbose = TRUE ; break ;
	case '?':
	    ERROR("Unrecognised option") ;
	}
    }
    
    /* Need to think whether this all works OK with indiv exclusions on top of 
       xinclude and yinclude, so */
    assert(indiv_exclude_file == NULL) ;
    n_inc = n_tot ;
    
/*     if(b <= a) { */
/* 	PRINT("a=%d >= b=%d!\n", a, b) ; */
/* 	exit(1) ; */
/*     } */
    assert(a >= 0 && b > a && c < d && d <= b) ;
    assert(b <= n_tot && d <= n_tot) ;
    
    if(c == a) {
	assert(d == b) ;
	symmetric = TRUE ;
    }
    else {
	symmetric = FALSE ;
	/* assert(c < d && d <= a) ;
	   Asserts that this code isn't being used to compute crossproducts for partially 
	   overlapping sets of individuals, and suggests that the lower diagonal of a larger
	   covariance matrix is being computed, in a blockwise fashion.
	   I don't really understand the above assertion any more, and it fails.
	*/
    }
    nx = b - a ;
    ny = d - c ;
    
    if(verbose) {
	PRINT("%s\tgenocov version %s\n", timestring(), version_string) ; fflush(stdout) ;
	PRINT("[a,b) = [%d,%d) -> nx = %d\n", a, b, nx) ;
	PRINT("[c,d) = [%d,%d) -> ny = %d\n", c, d, ny) ;
	PRINT("symmetric = %s\n", symmetric ? "TRUE" : "FALSE") ;
	PRINT("n_tot = %d\t", n_tot) ;
	PRINT("L_tot = %d\t", L_tot) ;
	PRINT("L_inc = %d\n", L_inc) ;
	PRINT("genotypes file = %s\n", genotypes_file) ;
	PRINT("data dir = %s\n", data_dir) ;
	PRINT("freqs file = %s\n", freqs_file) ;
	PRINT("indiv exclude file = %s\n", indiv_exclude_file) ;
	PRINT("SNP exclude file = %s\n", snp_exclude_file) ;
	PRINT("outdir = %s\n", outdir) ;
    }
    
    if(verbose) PRINT("setting indiv inclusion\n") ;
    indiv_include = make_include_vector(indiv_exclude_file, n_tot, n_inc) ;
    
    if(verbose) PRINT("setting SNP inclusion\n") ;
    snp_include = make_include_vector(snp_exclude_file, L_tot, L_inc) ;
    
    if(verbose) PRINT("reading freqs\n") ;
    freqs_inc = make_freqs_vector(freqs_file, L_tot, L_inc, snp_include) ;

    if(verbose) {
	PRINT("%s\t/* First %d entries of freqs: */\n", timestring(), MIN(10, L_inc)) ;
	for(i = 0 ; i < MIN(10, L_inc) ; i++) PRINT("%lf ", freqs_inc[i]) ;
	putchar('\n') ;
    }
    
    nmono = set_monomorphs_for_exclusion(snp_include, freqs_inc, L_tot, L_inc) ;
    
    if(nmono > 0) {
	if(verbose) PRINT("excluding %d monomorphs\n", nmono) ;
	L_inc -= nmono ;
	FREE(freqs_inc) ;
	freqs_inc = make_freqs_vector(freqs_file, L_tot, L_inc, snp_include) ;
    }
    
    /* 
       Set the names of the individual data files that will be needed,
       and construct individual inclusion vectors for the two subsets.
    */
    fnames = ALLOC(n_tot, char *) ;
    xinclude = ALLOC(n_tot, bool) ;
    nx = 0 ;
    if(!symmetric) {
	yinclude = ALLOC(n_tot, bool) ;
	ny = 0 ;
    }
    else {
	yinclude = NULL ;
	ny = nx ;
    }
    
    for(i_tot = 0 ; i_tot < n_tot ; i_tot++) {
	
	if(i_tot >= a && i_tot < b) {
	    fnames[i_tot] = ALLOC(100, char) ;
	    sprintf(fnames[i_tot], "%s/%05d", data_dir, i_tot + 1) ; // indexing from 1 in shell
	    xinclude[i_tot] = TRUE ;
	    nx++ ;
	}
	else xinclude[i_tot] = FALSE ;

	if(symmetric) continue ;
	else if(i_tot >= c && i_tot < d) {
	    fnames[i_tot] = ALLOC(100, char) ;
	    sprintf(fnames[i_tot], "%s/%05d", data_dir, i_tot + 1) ; // indexing from 1 in shell
	    yinclude[i_tot] = TRUE ;
	    ny++ ;
	}
	else yinclude[i_tot] = FALSE ;
    }
    
    if(verbose) PRINT("%s\treading genotypes for %d individuals at %d SNPs\n",
		      timestring(), nx, L_inc) ;
    
    if(genotypes_file != NULL) {
	// this results in x being n_tot x L_inc, column major
	x = ALLOC(L_inc * nx, double) ;
	fp = fopen(genotypes_file, "r") ;
	read_submatrix_genotypes_double_transpose(fp,
						  n_tot, xinclude, nx,
						  L_tot, snp_include, L_inc,
						  x) ;
	fclose(fp) ;
    }
    else x = read_genotypes(fnames, L_tot, L_inc, snp_include, n_tot, nx, xinclude, TRUE) ;
    
    if(verbose) {
	PRINT("%s\t/* First %d entries of genotypes x: */\n", timestring(), MIN(nx*L_inc,10)) ;
	for(i = 0 ; i < MIN(nx*L_inc,10) ; i++) PRINT("%lf ", x[i]) ;
	putchar('\n') ;
	PRINT("%s\tcentering and scaling %d x %d genotype matrix\n", timestring(), L_inc, nx) ;
    }
    centre_and_zero_missing_and_scale(x, L_inc, nx, freqs_inc) ;

    if(verbose) {
	PRINT("%s\t/* First %d entries of scaled genotypes x: */\n", timestring(), MIN(nx*L_inc,10)) ;
	for(i = 0 ; i < MIN(nx*L_inc,10) ; i++) PRINT("%lf ", x[i]) ;
	putchar('\n') ;
    }

    if(symmetric) {
	if(verbose) PRINT("%s\tforming symmetric crossproduct: t([%d,%d]) %%*%% [%d,%d]\n",
			  timestring(), L_inc, nx, L_inc, nx) ;
	cov = XTX(x, L_inc, nx, 1.0 / L_inc) ;
    }
    else {
	if(verbose) PRINT("%s\treading genotypes for %d individuals at %d SNPs\n", timestring(), nx, L_inc) ;
	if(genotypes_file != NULL) {
	    y = ALLOC(L_inc * ny, double) ;
	    fp = fopen(genotypes_file, "r") ;
	    read_submatrix_genotypes_double_transpose(fp,
						      n_tot, yinclude, ny,
						      L_tot, snp_include, L_inc,
						      y) ;
	    fclose(fp) ;
    	}
	else y = read_genotypes(fnames, L_tot, L_inc, snp_include, n_tot, ny, yinclude, TRUE) ;

	if(verbose) PRINT("%s\tcentering and scaling %d x %d genotype matrix\n",
			  timestring(), L_inc, ny) ;

	centre_and_zero_missing_and_scale(y, L_inc, ny, freqs_inc) ;
	
	if(verbose) PRINT("%s\tforming crossproduct: t([%d,%d]) %%*%% [%d,%d] \n",
			  timestring(), L_inc, nx, L_inc, ny) ;
	cov = matprod(x, y, TRUE, FALSE, L_inc, nx, L_inc, ny, 1.0 / L_inc) ;
	FREE(y) ;
	FREE(yinclude) ;
    }
    FREE(x) ;
    FREE(xinclude) ;
    FREE(snp_include) ;
    FREE(freqs_inc) ;
    
    if(verbose) PRINT("%s\twriting %d columns in dir %s/\n",
		      timestring(), symmetric ? nx : ny, outdir) ;

    for(j_tot = 0 ; j_tot < (symmetric ? nx : ny) ; j_tot++) {
	
	sprintf(outfile, "%s/%05d-%05d", outdir, j_tot + c + 1, a + 1) ;
	// a was xchunk; indexing from 1 in shell

	fp = fopen(outfile, "w") ; if(fp == NULL) ERROR("outfile open error\n") ;
	
	for(i_tot = (symmetric ? j_tot : 0) ; i_tot < nx ; i_tot++)
	    fprintf(fp, "%lf ", cov[i_tot + j_tot * nx]) ;

	if(b == n_tot) fputc('\n', fp) ;
	
	fclose(fp) ;
    }
    
    for(i_tot = 0 ; i_tot < n_tot ; i_tot++)
	if( (i_tot >= a && i_tot < b) || (i_tot >= c && i_tot < d) )
	    FREE(fnames[i_tot]) ;

    FREE(indiv_include) ;
    FREE(fnames) ;
    if(verbose) PRINT("%s\tdone\n", timestring()) ;

    return 0 ;
}
