#include <dan.h>
#include <unistd.h>
#include <assert.h>

int main(int argc, char *argv[]) {
    int c, bp, i, n, flipind ;
    double p[3] ;
    char snp_id[100], rs_id[100], alleles[2], *delim, *flipind_file ;
    FILE *f ;
    bool inconsistent ;

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
    
    f = fopen(flipind_file, "r") ;
    
    while( fscanf(stdin, "%s %s %d %c %c", snp_id, rs_id, &bp, alleles+0, alleles+1) == 5) {
	while((c = fgetc(f)) != '\n') {
	    assert(c != EOF) ;
	    flipind = c - '0' ; // so the flip indicator must be the last character on the line
	}
	printf("%s %s %d ", snp_id, rs_id, bp) ;
	if(flipind == 2) flipind = 0 ; // flipind == 2 indicates all missing data in one or other data set
	if(flipind != 0 && flipind != 1) {
	    assert(flipind == 3) ;
	    inconsistent = TRUE ;
	    flipind = 0 ;
	}
	else inconsistent = FALSE ;
	printf("%c %c ", alleles[flipind], alleles[1 - flipind]) ;

	for(i = 0 ; i < n ; i++) {
	    if( fscanf(stdin, "%lf %lf %lf", p+0, p+1, p+2) != 3)
		ERROR("fscanf error:") ;
	    if(inconsistent) printf("%d %d %d%s", 0, 0, 0, i+1 == n ? "" : " ") ;
	    else if(flipind) printf("%.2lf %.2lf %.2lf%s", p[2], p[1], p[0], i+1 == n ? "" : " ") ;
	    else printf("%.2lf %.2lf %.2lf%s", p[0], p[1], p[2], i+1 == n ? "" : " ") ;
	}
	putchar('\n') ;
    }
    
    fclose(f) ;
    return 0 ;
}
