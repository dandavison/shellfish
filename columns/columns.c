#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include <stdarg.h>

void print_stdout(char *fmt, ...) {
    va_list args;

    fflush(stdout);
    
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    va_end(args);
    
    fflush(stdout) ;
}
#define PRINT print_stdout
void ERROR (char *msg) {
    fputs(msg, stderr) ;
    fputc('\n', stderr) ;
    exit(2) ;
}

void WARN (char *msg) {
    fputs(msg, stderr) ;
    fputc('\n', stderr) ;
}

void usage() {
    ERROR("columns -f colnumfile < inputfile") ;
}

#define false 0
#define true 1
#define bool int

int main(int argc, char *argv[])
{
    
    size_t maxlinelength=0, maxidx=0;
    char *line=NULL, c;
    FILE *infile=stdin, *cfile=NULL;
    bool invert=false, strings=false;
    int *idx=NULL, *p, i, j, n, nin=-1, nout;

    while((c = getopt(argc, argv, "f:sv")) != -1) {
	switch(c) {
	case 'f':
	    cfile = fopen(optarg, "r") ;
	    break ;
	case 's':
	    strings = true ;
	    break ;
	case 'v':
	    invert = true;
	    break ;
	case '?':
	    usage() ;
	}
    }
    if(cfile == NULL)
	ERROR("Usage: columns -f indices [-s] [-v]");

    maxidx = 100;
    idx = (int *) malloc(maxidx * sizeof(int));
    nout = 0;
    while( fscanf(cfile, "%d", idx + nout) == 1 )
    {
	idx[nout]-- ; // indexing from one in shell
	nout++;
	if( nout == maxidx )
	{
	    maxidx *= 2;
	    idx = realloc(idx, maxidx * sizeof(int));
	}
    }

    for(i = 0 ; true ; ++i) {
	if(strings) {
	    n = process_line(idx, invert) ;
	}
	else 
	    n = getline(&line, &maxlinelength, infile) - 1 ; // newline char
	
	if(n < 0 ) break;
	if (i == 0) { // We're on the first line of input
	    nin = n;
	    // nout = invert ? nin - nout : nout;
	}
	else if(n != nin) {
	    fprintf(stderr, "Line 1 had length %d; line %d has length %d", nin, i+1, n);
	    ERROR("Input lines should be same length");
	}

	if(!strings) {
	    p = idx;
	    for(j = 0 ; j < nin ; j++) {
		if( (j != *p) == invert )
		    putchar(line[j]);
		if( j == *p ) ++p;
	    }
	}
	putchar('\n');
	if (i % 100 == 0)
	    fprintf(stderr, "%d\r", i);
    }
    return 0;
}


int process_line(int *idx, bool invert) {
    int rtn, i=0, j, *p ;
    char *buf=NULL ;
    size_t len ;
    
    for(p = idx, i = 0 ; true ; i++) {
	rtn = scanf("%ms", &buf) ;

	if(rtn == EOF) return -1 ;
	if(rtn != 1) ERROR("Failed to read a string, not EOF. This shouldn't happen") ;

	/* OK, we read another string */
	if( (i != *p) == invert )
	    printf("%s ", buf);
	if( i == *p ) ++p;
	free(buf) ;

	/* Now read following whitespace and check for end-of-line */
	rtn = scanf("%m[ \t\n\r\v\f]", &buf) ;
	len = strlen(buf) ;
	for(j = 0 ; j < len ; j++)
	    if(buf[j] == '\n' || buf[j] == '\r') {
		free(buf) ;
		return i ;
	    }
	free(buf) ;
    }
}
