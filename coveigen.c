#include <dan.h>
#include <io.h>
#include <unistd.h>
#include <assert.h>
#include <blas-lapack-dan.h>
#include <iogeno.h>

int main(int argc, char *argv[]) {
    int c, n, nc2, nvecs, i, j, k, v, rtn ;
    char *cov_file, *indiv_exclude_file, *outdir, outfile[500] ;
    double *s, *evals, *evecs ;
    FILE *fp ;
    bool verbose = FALSE;

    while((c = getopt(argc, argv, "n:V:g:I:o:v")) != -1) {
	switch(c) {
	case 'n':
	    n = atoi(optarg) ; break ;
	case 'V':
	    nvecs = atoi(optarg) ; break ;
	case 'g':
	    cov_file = optarg ; break ;
	case 'I':
	    indiv_exclude_file = optarg ; break ;
	case 'o':
	    outdir = optarg ; break ;
	case 'v':
	    verbose = TRUE ; break ;
	case '?':
	    ERROR("Unrecognised option") ;
	}
    }

    nc2 = n * (n - 1) / 2 ;
    s = ALLOC(n * n, double) ;
    fp = fopen(cov_file, "r") ;

    k = 0 ;
    for(j = 0 ; j < n ; j++) {
	if( verbose && !(j % 1000) ) PRINT("reading column %d of %d\n", j, n) ;
	for(i = j ; i < n ; i++) {
	    rtn = fscanf(fp, "%lf", s + i + n*j) ;
	    if(rtn != 1) {
		PRINT("return code = %d\n", rtn) ;
		ERROR("cov fscanf error\n") ;
	    }
	    k++ ;
	}
    }
    assert(k == nc2 + n) ;
    fclose(fp) ;

    evals = ALLOC(nvecs, double) ;
    evecs = ALLOC(n * nvecs, double) ;
    
    if(verbose) PRINT("%s\tdoing eigendecomposition\n", timestring()) ;
    eigen(s, n, nvecs, FALSE, evals, evecs) ;
    
    sprintf(outfile, "%s/evals", outdir) ;
    if(verbose) PRINT("%s\twriting evals to file %s\n", timestring(), outfile) ;
    fp = fopen(outfile, "w") ;
    write_evals(evals, fp, nvecs) ;
    fclose(fp) ;
    
    sprintf(outfile, "%s/evecs", outdir) ;
    if(verbose) PRINT("%s\twriting evecs to file %s\n", timestring(), outfile) ;
    fp = fopen(outfile, "w") ;
    write_evecs(evecs, fp, n, nvecs) ;
    fclose(fp) ;
    
    if(verbose) PRINT("%s\tdone!\n", timestring()) ;
    return 0 ;
}
