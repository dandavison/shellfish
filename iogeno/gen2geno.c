#include <dan.h>
#include <unistd.h>
#include <assert.h>

int main(int argc, char *argv[]) {
    int L, l, n, i, c, bp, geno ;
    double thresh, p[3] ;
    char snp_id[100], rs_id[100], alleles[2], *delim ;
    
    delim="" ;
    while((c = getopt(argc, argv, "n:t:d:")) != -1) {
	switch(c) {
	case 't':
	    thresh = atof(optarg) ; break ;
	case 'n':
	    n = atoi(optarg) ; break ;
	case 'd':
	    delim = optarg ; break ;
	case '?':
	    ERROR("unrecognised option\n") ;
	}
    }
    
    assert(thresh > 0.5) ;

    while( fscanf(stdin, "%s %s %d %c %c", snp_id, rs_id, &bp, alleles+0, alleles+1) == 5) {
	for(i = 0 ; i < n ; i++) {
	    if( fscanf(stdin, "%lf %lf %lf", p+0, p+1, p+2) != 3)
		ERROR("fscanf error:") ;
	    geno = p[0] >= thresh ? 0 : p[1] >= thresh ? 1 : p[2] >= thresh ? 2 : MISSING ;
	    if(i > 0) printf("%s", delim) ;
	    putchar('0' + geno) ;
	}
	putchar('\n') ;
    }
    
    return 0 ;
}
