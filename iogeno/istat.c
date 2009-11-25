#include <dan.h>
#include <assert.h>
#include <unistd.h>

void usage() {
    ERROR("istat -n n -m -p") ;
}

void update_missing_rate_chiamo(FILE *f, int n, double thresh, double *x) {
    int i ;
    double p[3] ;
    
    for(i = 0 ; i < n ; i++) {
	if( fscanf(stdin, "%lf %lf %lf", p+0, p+1, p+2) != 3)
	    ERROR("fscanf error:") ;
	x[i] += 1 - p[0] - p[1] - p[2] ;
    }
}

int main(int argc, char *argv[]) {
    unsigned char *snp ;
    int i, n, c, bp, L ;
    char stat, snp_id[100], rs_id[100], alleles[2] ;
    bool chiamo ;
    double *x, thresh ;
    
    n = -1 ;
    stat = '?' ;
    chiamo = FALSE ;
    
    while((c = getopt(argc, argv, "mn:pt:")) != -1) {
	switch(c) {
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
    if(n < 0) ERROR("arguments: n\n") ;
    if(stat == '?') usage() ;
    
    x = CALLOC(n, double) ;
    
    L = 0 ;
    if(chiamo)
	while( fscanf(stdin, "%s %s %d %c %c", snp_id, rs_id, &bp, alleles+0, alleles+1) == 5 ) {
	    switch(stat) {
	    case 'm':
		update_missing_rate_chiamo(stdin, n, thresh, x) ; break ;
	    }
	    L++ ;
	}
    else {
	snp = ALLOC(n, unsigned char) ;
	while(fread(snp, 1, n, stdin) == n) {
	    for(i = 0 ; i < n ; ++i)
		if( snp[i] - '0' == MISSING ) x[i]++ ;
	    L++;
	}
	FREE(snp) ;
    }
    
    for(i = 0 ; i < n ; i++)
	printf("%lf\n", x[i] / L) ;
    
    return 0 ;
}
