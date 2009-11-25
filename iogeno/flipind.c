#include <dan.h>
#include <unistd.h>
#include <assert.h>

#define NA '-'
#define ALL_MISSING 2
#define INCONSISTENT 3

int validate(char *s1, char *s2, int snp_num) {
    char alleles[5] ;
    bool ok1, ok2 ;
    int i, j, ans = -1 ;

    alleles[0] = 'A' ;
    alleles[1] = 'C' ;
    alleles[2] = 'G' ;
    alleles[3] = 'T' ;
    alleles[4] = '-' ;

    for(i = 0 ; i < 2 ; i++) {
	ok1 = FALSE ;
	ok2 = FALSE ;
	for(j = 0 ; j <= 4 ; j++) {
	    if(s1[i] == alleles[j]) ok1 = TRUE ;
	    if(s2[i] == alleles[j]) ok2 = TRUE ;
	}
	assert(ok1) ;
	assert(ok2) ;
    }

    if(s1[0] == 'A' && s1[1] == 'T' ||
       s1[0] == 'T' && s1[1] == 'A' ||
       s2[0] == 'A' && s2[1] == 'T' ||
       s2[0] == 'T' && s2[1] == 'A' ||

       s1[0] == 'C' && s1[1] == 'G' ||
       s1[0] == 'G' && s1[1] == 'C' ||
       s2[0] == 'C' && s2[1] == 'G' ||
       s2[0] == 'G' && s2[1] == 'C') {
	/* EPRINT("SNP %d: no A/T or G/C SNPs allowed: SNP1 = %c/%c and SNP2 = %c/%c\n", */
	/*        snp_num, s1[0], s1[1], s2[0], s2[1]) ; */
	ans = INCONSISTENT ;
    }
    return ans ;
} 

char flip(char allele) {
    if(allele == 'A') return 'T' ;
    if(allele == 'C') return 'G' ;
    if(allele == 'G') return 'C' ;
    if(allele == 'T') return 'A' ;
    ERROR("flip error: allele = %c", allele) ;
}


int main(int argc, char *argv[]) {
    char s1[2], s2[2], delim='\t' ;
    int ans ;
    int i = 0, bad = 0 ;

    //    while(fscanf(stdin, "%c %c %c %c", s1+0, s1+1, s2+0, s2+1) == 4) {
    while(fscanf(stdin, "%c\t%c\t%c\t%c", s1+0, s1+1, s2+0, s2+1) == 4) {
	i++ ;
	assert(getchar() == '\n') ;
	
	// PRINT("%c%c|%c%c ", s1[0], s1[1], s2[0], s2[1]) ;
	
	if( validate(s1, s2, i) == INCONSISTENT ) {
	    bad++ ;
	    PRINT("%d\n", INCONSISTENT) ;
	    continue ;
	}
	if(s1[0] == NA) { // all missing
	    assert(s1[1] == NA) ;
	    PRINT("%d\n", ALL_MISSING) ;
	    continue ;
	}
	if(s2[0] == NA) { // all missing
	    assert(s2[1] == NA) ;
	    PRINT("%d\n", ALL_MISSING) ;
	    continue ;
	}
	if(s1[1] == NA && s2[1] == NA) { // both are monomorphic
	    if( !(s1[0] == s2[0] || flip(s1[0]) == s2[0]) ) {
		EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
		EPRINT("SNP %d: both are monomorphic, but alleles are inconsistent", i) ;
		// exit(2) ;
		bad++ ;
		PRINT("%d\n", INCONSISTENT) ;
		
		continue ;
	    }
	    PRINT("%d\n", s1[0] != s2[0]) ;
	    continue ;
	}
	if(s1[1] == NA) { // first is monomorphic, second isn't
	    if(s1[0] == s2[0] || flip(s1[0]) == s2[0]) ans = 0 ;
	    else if(s1[0] == s2[1] || flip(s1[0]) == s2[1]) ans = 1 ;
	    else {
		EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
		EPRINT("SNP %d: first is monomorphic, but alleles are inconsistent", i) ;
		// exit(2) ;
		bad++ ;
		ans = INCONSISTENT ;
	    }
	    PRINT("%d\n", ans) ;
	    continue ;
	}
	if(s2[1] == NA) { // second is monomorphic, first isn't
	    if(s2[0] == s1[0] || flip(s2[0]) == s1[0]) ans = 0 ;
	    else if(s2[0] == s1[1] || flip(s2[0]) == s1[1]) ans = 1 ;
	    else {
		EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
		EPRINT("SNP %d: second is monomorphic, but alleles are inconsistent", i) ;
		// exit(2) ;
		bad++ ;
		ans = INCONSISTENT ;
	    }
	    PRINT("%d\n", ans) ;
	    continue ;
	}
	
	/* If we get here, both SNPs have two alleles. */
	assert(s1[0] != s1[1] && s2[0] != s2[1]) ;
	
	if(s1[0] == s2[0]) { // same strand, same reference allele
	    if(s1[1] != s2[1]) {
		EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
		EPRINT("SNP %d: since s1[0] == s2[0], expecting s1[1] == s2[1]\n", i) ;
		// exit(2) ;
		bad++ ;
		ans = INCONSISTENT ;
	    }
	    else ans = 0 ;
	    PRINT("%d\n", ans) ;
	    continue ;
	}
	/* So the first alleles differ: reference allele differs, or strand differs, or both */
	
	if(s1[1] == s2[0]) { // same strand, different reference allele
	    // (we must be on the same strand, because A/T and C/G SNPs are not allowed)
	    if(s1[0] != s2[1]) {
		EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
		EPRINT("SNP %d: since s1[1] == s2[0], expecting s1[0] == s2[1]\n", i) ;
		// exit(2) ;
		bad++ ;
		ans = INCONSISTENT ;
	    }
	    else ans = 1 ;
	    PRINT("%d\n", ans) ;
	    continue ;
	}
	if(flip(s1[0]) == s2[0]) { // different strand, but same reference allele
	    if(flip(s1[1]) != s2[1]) {
		EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
		EPRINT("SNP %d: since the 1st agree after flipping, expecting 2nd alleles to do so\n", i) ;
		// exit(2) ;
		bad++ ;
		ans = INCONSISTENT ;
	    }
	    else ans = 0 ;
	    PRINT("%d\n", ans) ;
	    continue ;
	}
	if(flip(s1[0]) == s2[1]) { // different strand, different reference allele
	    if(flip(s1[1]) != s2[0]) {
		EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
		EPRINT("SNP %d: since 1st and 2nd agree after flipping, expecting 2nd and 1st to do so\n", i) ;
		// exit(2) ;
		bad++ ;
		ans = INCONSISTENT ;
	    }
	    else ans = 1 ;
	    PRINT("%d\n", ans) ;
	    continue ;
	}

	EPRINT("%c %c | %c %c\t", s1[0], s1[1], s2[0], s2[1]) ;
	EPRINT("SNP %d: these allele pairs seem inconsistent\n", i) ;
	// exit(2) ;
	bad++ ;
	PRINT("%d\n", INCONSISTENT) ;
    }
    if( bad > 0 ) EPRINT("%d inconsistent SNPs out of %d\n", bad, i) ;
    return 0 ;
}
