#ifndef _HAVE_DAN_H
#define _HAVE_DAN_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#define _HAVE_BLAS_LAPACK 1

// #define MATHLIB_STANDALONE TRUE
// #include <Rmath.h>

void print_stdout(char *fmt, ...) ;
void print_stderr(char *fmt, ...) ;
void eprintf(char *fmt, ...) ;
void weprintf(char *fmt, ...) ;
#define PRINT print_stdout
#define EPRINT print_stderr
#define ERROR eprintf
#define WARN weprintf

void *malloc_dan(size_t size) ;
void *calloc_dan(size_t nmemb, size_t size) ;

#define ALLOC(howmany, type) (type *) malloc_dan((howmany) * sizeof(type))
#define CALLOC(howmany, type) (type *) calloc_dan((howmany), sizeof(type))
#define FREE(ptr) free(ptr)


char *timestring() ;

#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))

typedef int bool ;
#define TRUE 1
#define FALSE 0

#define true 1
#define false 0

#define MISSING 9


#define INDIV_MAJOR_ORDER TRUE /* genotypes stored as if {SNPs x indivs} array, in column major order,
				  i.e. all SNPs for first indiv, then all SNPs for second, ... */


#endif // _HAVE_DAN_H
