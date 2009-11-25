/* 
   Template functions to avoid duplicating io code for reading different types.

   Immediately prior to #include-ing these functions, you need to
   #define TYPE as double, int, etc, and also #define each of the
   function names in this file accordingly. E.g.

   #define TYPE double
   #define read_matrix read_matrix_double
   #define write_matrix write_matrix_double

   #include "io-templates.c"
   
   #undef TYPE
   #undef read_matrix
   #undef write_matrix

   This is done in io.c, so if io.c is included in the build then you shouldn't 
   need to include this file.

*/

void read_submatrix(FILE *f, int majordim, bool *majorinc, int minordim, bool *minorinc,
		    char *format, TYPE *dest) {
    
    /*
      f is pointing at the start of formatted text corresponding to majordim x minordim elements of a matrix.
      e.g. if the matrix is in column-major order, then majordim is the number of rows, and minordim the
      number of columns.
      
      The majorinc and minorinc have lengths majordim and minordim respectively with element[i] indicating
      whether to include the corresponding row/column.
    
      Thus although majordim is necessarily the actual major dimension of the matrix on file, minordim may be
      less than the minor dimension (and that actual minor dimension need not be known by this function).
      
    */
    
    int i, j, rtn, nout=0 ;
    bool output_line, output_col ;
    TYPE temp;

    for(j = 0 ; j < minordim ; j++) {
	output_line = (minorinc == NULL || minorinc[j]) ;
	for(i = 0 ; i < majordim ; i++) {
	    rtn = fscanf(f, format, &temp) ;
	    if( rtn != 1 ) ERROR("read_submatrix_TYPE: fscanf read %d elements error:", rtn) ;
	    output_col = (majorinc == NULL || majorinc[i]) ;
	    if(output_line && output_col)
		dest[nout++] = temp ;
	}
    }
}

void read_submatrix_transpose(FILE *f,
			      int majordim, bool *majorinc, int majorincn,
			      int minordim, bool *minorinc, int minorincn,
			      char *format, TYPE *dest) {
    
    /*
      The same as read_submatrix except that output is transposed,
    */

    int i, j, iout, jout, rtn ;
    bool output_column ;
    TYPE temp ;

    for(j = jout = 0 ; j < minordim ; j++) {
	output_column = ( minorinc == NULL || minorinc[j] ) ;
	for(i = iout = 0 ; i < majordim ; i++) {
	    rtn = fscanf(f, format, &temp) ;
	    if( rtn != 1 ) ERROR("read_submatrix_transpose_TYPE: fscanf read %d elements error:", rtn) ;
	    if( output_column && (majorinc == NULL || majorinc[i]) )
		dest[iout++ * minorincn + jout] = temp ; // (row jout, col iout)
	}
	if( output_column ) jout++ ;
    }
}

void write_matrix(TYPE *p, FILE *f, int majordim, int minordim, char *format) {
    int i, j ;
    
    for(j = 0 ; j < minordim ; j++) {
	for(i = 0 ; i < majordim ; i++)
	    fprintf(f, format, *p++) ;
	fprintf(f, "\n") ;
    }
}

TYPE *read_matrix(char *fname, int dim1, int dim2, char *format) {
    FILE *f ;
    TYPE *p, *p0 ;
    int i, j, rtn ;
    char c ;

    f = fopen(fname, "r") ;
    if(f == NULL) ERROR("read_matrix_TYPE: fopen error\n") ;
    
    p = p0 = ALLOC(dim1 * dim2, TYPE) ;
    
    for(i = 0 ; i < dim1 ; i++) {
	for(j = 0 ; j < dim2 ; j++)
	    if( (rtn = fscanf(f, format, p++)) != 1 ) {
		WARN("return code = %d\n", rtn) ;
		ERROR("read_matrix_TYPE fscanf error\n") ;
	    }
	while( (c = fgetc(f)) == ' ') ;
	if(c != '\n') ERROR("expecting newline (perhaps after spaces) but read '%c'", c) ;
    }
    assert(fgetc(f) == EOF) ;
    
    fclose(f) ;
    return p0 ;
}
