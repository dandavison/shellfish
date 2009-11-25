#include <dan.h>

void *malloc_dan(size_t size) {
    void *p = malloc(size) ;
    if(p == NULL) ERROR("malloc failed:") ;
    return p ;
}

void *calloc_dan(size_t nmemb, size_t size) {
    void *p = calloc(nmemb, size) ;
    if(p == NULL) ERROR("calloc failed:") ;
    return p ;
}

char *timestring() {
    static char buff[100] ;
    time_t t ;
    time(&t) ;
    strftime(buff, 100, "%a %d %b %H:%M:%S", localtime(&t)) ;
    return buff ;
}


/* 
   Haven't been bothered to work out how to write a single function that takes a file pointer and a variable argument
   list and then redefine the following two in terms of calls to that function.
 */

void print_stdout(char *fmt, ...) {
    va_list args;

    fflush(stdout);
    
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    va_end(args);
    
    fflush(stdout) ;
}

void print_stderr(char *fmt, ...) {
    va_list args;

    fflush(stderr);
    
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    
    fflush(stderr) ;
}

// Following functions derive from

/* Copyright (C) 1999 Lucent Technologies */
/* Excerpted from 'The Practice of Programming' */
/* by Brian W. Kernighan and Rob Pike */

//  Modification of eprintf
void eprintf(char *fmt, ...)
{
    va_list args;

    fflush(stdout);
    
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    if (fmt[0] != '\0' && fmt[strlen(fmt)-1] == ':')
	fprintf(stderr, " %s", strerror(errno));
    fprintf(stderr, "\n");
    exit(2); /* conventional value for failed execution */
}

/* weprintf: print warning message */
void weprintf(char *fmt, ...)
{
	va_list args;

	fflush(stdout);
	fprintf(stderr, "warning: ");

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	if (fmt[0] != '\0' && fmt[strlen(fmt)-1] == ':')
		fprintf(stderr, " %s\n", strerror(errno));
	else
		fprintf(stderr, "\n");
}
