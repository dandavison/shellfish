#include <dan.h>
#include <assert.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    unsigned char *snp ;
    char * flipind_file ;
    FILE *f ;
    int n, c, l, i, geno, flipind ;
    
    n = -1 ;
    while((c = getopt(argc, argv, "i:n:")) != -1) {
	switch(c) {
	case 'i':
	    flipind_file = optarg ; break ;
	case 'n':
	    n = atoi(optarg) ; break ;
	case '?':
	    ERROR("unrecognised option\n") ;
	}
    }
    if(n < 0) ERROR("arguments: n\n") ;
    
    snp = ALLOC(n, unsigned char) ;
    f = fopen(flipind_file, "r") ;
    
    l=1 ;
    while(fread(snp, 1, n, stdin) == n) {
	while((c = fgetc(f)) != '\n') {
	    assert(c != EOF) ;
	    flipind = c - '0' ; // so the flip indicator must be the last character on the line
	}
	if(flipind != 0 && flipind != 1) {
	    // fprintf(stderr, "flipind at locus %d is %d (should be 0/1)\n", l, flipind) ;
	    // exit(2) ;
	    geno = MISSING ;
	}
	for(i = 0 ; i < n ; i++) {
	    geno = snp[i] - '0' ;
	    assert(geno >= 0 && geno <= 2 || geno == MISSING) ; 
	    if(flipind == 1 && geno != MISSING)
		geno = 2 - geno ;
	    putchar('0' + geno) ;
	}
	assert(getchar() == '\n') ;
	putchar('\n') ;
	l++ ;
    }
    FREE(snp) ;
    return 0 ;
}
