#include <dan.h>
#include <assert.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    unsigned char *snp ;
    int n, c, g, l, i ;
    bool probs_only ;
    
    n = -1 ;
    probs_only = FALSE ;
    while((c = getopt(argc, argv, "n:p")) != -1) {
	switch(c) {
	case 'n':
	    n = atoi(optarg) ; break ;
	case 'p':
	    probs_only = TRUE ; break ;
	case '?':
	    ERROR("unrecognised option\n") ;
	}
    }
    if(n < 0) ERROR("arguments: n\n") ;
    
    snp = ALLOC(n, unsigned char) ;
    
    l=1 ;
    while(fread(snp, 1, n, stdin) == n) {
	if(!probs_only) printf("snp_%06d rs%06d 00000000 - - ", l, l) ;
	for(i = 0 ; i < n ; i++) {
	    g = snp[i] - '0' ;
	    printf("%d %d %d%s", g == 0, g == 1, g == 2, i+1 == n ? "" : " ") ;
	}
	assert(getchar() == '\n') ;
	putchar('\n') ;
	l++ ;
    }
    FREE(snp) ;
    return 0 ;
}
