#include <dan.h>
#include <assert.h>
#include <unistd.h>

double compute_freq_geno(unsigned char *snp, int n) {
    int g, i, x, n_obs ;

    n_obs = x = 0 ;
    for(i = 0 ; i < n ; i++) {
	g = snp[i] - '0' ;
	if(g != MISSING) {
	    x += g ;
	    n_obs++ ;
	}
    }
    return x / (2.0 * n_obs) ;
}

double compute_missing_rate_geno(unsigned char *snp, int n) {
    int i, count ;
    
    count  = 0 ;
    for(i = 0 ; i < n ; i++)
	if(snp[i] - '0' == MISSING) count++ ;
    
    return count / (double) n ;
}


double compute_freq_chiamo(FILE *f, int n, double thresh) {
    int i, rtn ;
    double p[3], pnum, ptot ;
    
    pnum = ptot = 0.0 ;
    for(i = 0 ; i < n ; i++) {
	rtn = fscanf(stdin, "%lf %lf %lf", p+0, p+1, p+2) ;
	if(rtn != 3) ERROR("fscanf failed to read requested number of items") ;
	pnum += 0.5 * p[1]  +  p[2] ;
	ptot += p[0] + p[1] + p[2] ;
    }

    return pnum / ptot ;
}

double compute_missing_rate_chiamo(FILE *f, int n, double thresh) {
    int i, rtn ;
    double p[3], pnum ;
    
    pnum = 0.0 ;
    for(i = 0 ; i < n ; i++) {
	rtn = fscanf(stdin, "%lf %lf %lf", p+0, p+1, p+2) ;
	if(rtn != 3) ERROR("fscanf failed to read requested number of items") ;
	pnum += 1 - p[0] - p[1] - p[2] ;
    }

    return pnum / n ;
}

int main(int argc, char *argv[]) {
    unsigned char *snp ;
    int n, c, bp ;
    char stat, snp_id[100], rs_id[100], alleles[2] ;
    bool chiamo, keep_legend ;
    double x, thresh = -1 ;
    
    n = -1 ;
    stat = 'f' ;
    chiamo = FALSE ;
    keep_legend = FALSE ;

    while((c = getopt(argc, argv, "kmn:pt:")) != -1) {
	switch(c) {
	case 'k':
	    keep_legend = TRUE ; break ;
	case 'm':
	    stat = 'm' ; break ;
	case 'n':
	    n = atoi(optarg) ; break ;
	case 'p':
	    chiamo = TRUE ; break ;
	case 't':
	    ERROR("Nothing implemented yet for chiamo thresh option -t") ;
	    thresh = atof(optarg) ;
	    break ;
	case '?':
	    ERROR("unrecognised option\n") ;
	}
    }
    if(n < 0) ERROR("usage: sstat -n num_indivs [-mpk]\n") ;
    
    if(chiamo)
	while( fscanf(stdin, "%s %s %d %c %c", snp_id, rs_id, &bp, alleles+0, alleles+1) == 5 ) {
	    if(keep_legend) PRINT("%s %s %d %c %c ", snp_id, rs_id, bp, alleles[0], alleles[1]) ;
	    switch(stat) {
	    case 'f':
		x = compute_freq_chiamo(stdin, n, thresh) ; break ;
	    case 'm':
		x = compute_missing_rate_chiamo(stdin, n, thresh) ; break ;
	    }
	    printf("%lf\n", x) ;
 	}
    else {
	snp = ALLOC(n, unsigned char) ;
	while(fread(snp, 1, n, stdin) == n) {
	    switch(stat) {
	    case 'f':
		x = compute_freq_geno(snp, n) ; break ;
	    case 'm':
		x = compute_missing_rate_geno(snp, n) ; break ;
	    }
	    printf("%lf\n", x) ;
	    assert(getchar() == '\n') ;
	}
	FREE(snp) ;
    }
    
    return 0 ;
}
