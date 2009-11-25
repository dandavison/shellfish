#include <dan.h>
#include <io.h>
#include <assert.h>

/* NB If you still think there's any point in C++ after seeing this
   file then there's something wrong with you. :) */

#define TYPE double
#define read_submatrix read_submatrix_double
#define read_submatrix_transpose read_submatrix_transpose_double
#define read_matrix read_matrix_double
#define write_matrix write_matrix_double
#include <io-templates.c>

#undef TYPE
#undef read_submatrix
#undef read_submatrix_transpose
#undef read_matrix
#undef write_matrix

#define TYPE int
#define read_submatrix read_submatrix_int
#define read_submatrix_transpose read_submatrix_transpose_int
#define read_matrix read_matrix_int
#define write_matrix write_matrix_int
#include <io-templates.c>

#undef TYPE
#undef read_submatrix
#undef read_submatrix_transpose
#undef read_matrix
#undef write_matrix

#define TYPE unsigned char
#define read_submatrix read_submatrix_char
#define read_submatrix_transpose read_submatrix_transpose_char
#define read_matrix read_matrix_char
#define write_matrix write_matrix_char
#include <io-templates.c>

#undef TYPE
#undef read_submatrix
#undef read_submatrix_transpose
#undef read_matrix
#undef write_matrix
