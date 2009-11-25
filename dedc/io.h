#ifndef _HAVE_IO_H
#define _HAVE_IO_H

void write_matrix_double(double *p, FILE *f, int dim1, int dim2, char *format) ;
double *read_matrix_double(char *fname, int dim1, int dim2, char *format) ;
void read_submatrix_double(FILE *f, int majordim, bool *majorinc, int minordim, bool *minorinc, char *format, double *dest) ;
//void read_submatrix_transpose_double(FILE *f, int majordim, bool *majorinc, int minordim, bool *minorinc, char *format, double *dest) ;

void write_matrix_int(int *p, FILE *f, int dim1, int dim2, char *format) ;
int *read_matrix_int(char *fname, int dim1, int dim2, char *format) ;

#endif // _HAVE_IO_H
