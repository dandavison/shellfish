#include <math.h>
#include <assert.h>

#include <dan.h>
#include <io.h>


void read_chiamo(FILE *f, int n, int L_tot, bool *snp_include, double *dest) {
    int i, l, bp, rtn ;
    char snp_id[100], rs_id[100], alleles[2] ;
    
    for(l = 0 ; l < L_tot ; l++) {
	rtn = fscanf(f, "%s %s %d %c %c", snp_id, rs_id, &bp, alleles+0, alleles+1) ;
	assert(rtn == 5) ;
	for(i = 0 ; i < n ; i++) {
	    if( fscanf(f, "%lf %lf %lf", dest+0, dest+1, dest+2) != 3 )
		ERROR("read_chiamo: fscanf error:") ;
	    if(snp_include[l]) dest += 3 ;
	}
    }
}


/* 
   Files containing genotypes have one line per SNP, 
   and no spaces.
*/

void read_submatrix_genotypes_double(FILE *f, int majordim, bool *majorinc, int minordim,
				     bool *minorinc, double *dest) {

    /*
      f is pointing at the start of text corresponding to majordim x minordim elements of a matrix,
      all of which are integers in the range 0-9, and between which there are no spaces.
      After each majordim entries there is a newline character.
      
      e.g. if the matrix is in column-major order, then majordim is the number of rows,
      and minordim the number of columns.
      
      The majorinc and minorinc have lengths majordim and minordim respectively with element[i]
      indicating whether to include the corresponding row/column.
      
      Thus although majordim is necessarily the actual major dimension of the matrix on file,
      minordim may be less than the minor dimension (and that actual minor dimension need not be
      known by this function).
    */

    int i, j;
    char c ;
    unsigned char *buff = ALLOC(majordim, unsigned char) ;
    bool seekable ;

    seekable = ( fseek(f, 0, SEEK_CUR) == 0 ) ; // if(!seekable) stdin presumably a pipe
    
    for(j = 0 ; j < minordim ; j++) {
	if(minorinc == NULL || minorinc[j]) {
	    if( fread(buff, 1, majordim, f) != majordim )
		ERROR("read_submatrix_genotypes_double: fread error:") ;
	    for(i = 0 ; i < majordim ; i++)
		if(majorinc == NULL || majorinc[i])
		    *dest++ = buff[i] - '0' ;
	}
	else {
	    if(seekable) {
		if( fseek(f, majordim, SEEK_CUR) != 0 )
		    ERROR("read_submatrix_genotypes_double: fseek error:") ;
	    }
	    else if( fread(buff, 1, majordim, f) != majordim )
		ERROR("read_submatrix_genotypes_double: fread error:") ;
	}
	if((c = fgetc(f)) != '\n') // NB not just a check; file pointer advance is necessary
	    ERROR("expecting newline, got %c (j = %d)", c, j) ; 
    }
    
    FREE(buff) ;
}

void read_submatrix_genotypes_double_transpose(FILE *f,
					       int majordim, bool *majorinc, int majorincn,
					       int minordim, bool *minorinc, int minorincn,
					       double *dest) {
    
    /*
      The same as read.submatrix_genotypes_double except that output is transposed,
    */

    int i, j, iout, jout ;
    char c ;
    unsigned char *buff = ALLOC(majordim, unsigned char) ;
    bool seekable ;

    seekable = ( fseek(f, 0, SEEK_CUR) == 0 ) ; // if(!seekable) stdin presumably a pipe
    
    for(j = jout = 0 ; j < minordim ; j++) {
	if(minorinc == NULL || minorinc[j]) {
	    if( fread(buff, 1, majordim, f) != majordim )
		ERROR("read_submatrix_genotypes_double: fread error:") ;
	    for(i = iout = 0 ; i < majordim ; i++)
		if(majorinc == NULL || majorinc[i]) {
		    dest[iout++ * minorincn + jout] = buff[i] - '0' ; // (row jout, col iout)
		    // non-transpose verion was *dest++ = buff[i] - '0' ;
		    // which is equiv to dest[jout * majornincn + iout] (row iout, col jout)
		}
	    jout++ ;
	}
	else {
	    if(seekable) {
		if( fseek(f, majordim, SEEK_CUR) != 0 )
		    ERROR("read_submatrix_genotypes_double: fseek error:") ;
	    }
	    else if( fread(buff, 1, majordim, f) != majordim )
		ERROR("read_submatrix_genotypes_double: fread error:") ;
	}
	if((c = fgetc(f)) != '\n') // NB not just a check; file pointer advance is necessary
	    ERROR("expecting newline, got %c (j = %d)", c, j) ; 
    }
    
    FREE(buff) ;
}

void centre_and_zero_missing_and_scale(double *g, int L, int n, double *p) {
    int l, i ;
    
    /* 
       Patterson et al procedure prior to forming crossproduct. 
       There may be a fancier way to do this with BLAS, but I'll
       leave it for now as I don't think it takes long to do like this.
    */
    
    assert(INDIV_MAJOR_ORDER) ;
    
    for(i = 0 ; i < n ; i++)
	for(l = 0 ; l < L ; l++) {
	    if(g[l + i*L] != MISSING) {
		g[l + i*L] -= 2.0 * p[l] ;
		g[l + i*L] /= sqrt(2 * p[l] * (1 - p[l])) ;
	    }
	    else g[l + i*L] = 0.0 ;
	}
}

unsigned char *read_genotypes_snp_range_all_indivs(char *file, int n, int snp_start, int snp_stop) {
    FILE *f = fopen(file, "r") ;
    unsigned char *x, *x0 ;
    int l, L = snp_stop - snp_start ;
    size_t nread ;

    x = x0 = ALLOC(n * L, unsigned char) ;
    
    // go to beginning of first line of data in the requested SNP range (n+1 because newlines)
    fseek(f, snp_start * (n+1), SEEK_SET) ;

    for(l = 0 ; l < L ; l++) {
	if(fread(x, 1, n, f) != n)
	    ERROR("fread failed to read n items");
	
	// can't comment out next line, obviously, as file pointer advance is necessary
	if(fgetc(f) != '\n')
	    ERROR("error reading genotypes: are you sure n was specified correctly?") ;
	x += n ;
    }
    
    fclose(f) ;

    return x0 ;
}


int *read_individual(int i, FILE *f, int n, int L_tot, int L_inc, bool *include) {
    /* Read column i from the genotype file */
    int c, l_tot, l_inc ;
    int *x = ALLOC(L_inc, int) ;
    
    fseek(f, i, SEEK_SET) ; // jump to requested column in first line
    l_inc = 0 ;
    for(l_tot = 0 ; l_tot < L_tot ; l_tot++) {
	c = fgetc(f) ;
	if(include[l_tot]) x[l_inc++] = c - '0' ;
	fseek(f, n, SEEK_CUR) ;
    }
    return x ;
}


/*
  Following functions deal with files containing genotypes for single individuals
*/

void read_individual_alone(char *fname, int L_tot, int L_inc, bool *include, double *dest) {
    /* Read individual from its own genotype file */
    int l_tot, l_inc ;
    FILE *f ;
    unsigned char c ;
    
    f = fopen(fname, "r") ;
    if(f == NULL) {
	PRINT("Trying to read %s\n", fname) ;
	ERROR("read_individual_alone: fopen error\n") ;
    }
    
    l_inc = 0 ;
    for(l_tot = 0 ; l_tot < L_tot ; l_tot++) {
	c = fgetc(f) ;
	if(include[l_tot]) {
	    dest[l_inc++] = c - '0' ;
	}
    }
    assert(fgetc(f) == '\n') ;
    fclose(f) ;
}

double *read_genotypes(char **fnames,
		       int L_tot, int L_inc, bool *snp_include,
		       int n_tot, int n_inc, bool *indiv_include, bool indiv_major_order) {
    int i_inc, i_tot ;
    if(!indiv_major_order) ERROR("read_genotypes only implemented for individual-major order\n") ;
    
    double *g = ALLOC(n_inc * L_inc, double) ;
    
    i_inc = 0 ;
    for(i_tot = 0 ; i_tot < n_tot ; i_tot++)
      if(indiv_include[i_tot]) {
	read_individual_alone(fnames[i_tot], L_tot, L_inc, snp_include, g + i_inc * L_inc) ;
	i_inc++ ;
      }
    
    return g ;
}

int *read_whole_individual_alone(FILE *f, int L) {
    /* Read column individual genotype file */
    int l ;
    int *x = ALLOC(L, int) ;
    if(fread(x, 1, L, SEEK_SET) != L)
	ERROR("fread failed to read requested number of items") ;
    
    for(l = 0 ; l < L ; l++)
	x[l] = x[l] - '0' ;
    
    return x ;
}

unsigned char *read_whole_individual_char(int i, FILE *f, int n, int L) {
    /* Read column i from the genotype file */
    int  l ;
    unsigned char *x = ALLOC(L, unsigned char) ;
    
    fseek(f, i, SEEK_SET) ; // jump to requested column in first line
    
    for(l = 0 ; l < L ; l++) {
	x[l] = fgetc(f) ;
	fseek(f, n, SEEK_CUR) ;
    }
    return x ;
}


/******  END SINGLE-INDIVIDUAL FUNCTIONS  ********/


bool *make_include_vector(char *exclude_file, int L_tot, int L_inc) {
    int L_exc = L_tot - L_inc, l_tot, l_exc ;
    int *exclude ;
    bool *include ;
    
    if(exclude_file != NULL)
	exclude = read_matrix_int(exclude_file, L_exc, 1, "%d") ;
    else {
	assert(L_exc == 0) ;
	exclude = NULL ;
    }
    
    include = ALLOC(L_tot, bool) ;
    
    /* Check the exclude indices are strictly increasing */
    for(l_exc = 1 ; l_exc < L_exc ; l_exc++)
	assert(exclude[l_exc] > exclude[l_exc-1]) ;
    
    l_exc = 0 ;
    for(l_tot = 0 ; l_tot < L_tot ; l_tot++) {
	if(L_exc > 0 && l_tot == (exclude[l_exc] - 1)) { // indices in exclude files start at 1
	    include[l_tot] = FALSE ;
	    l_exc++ ;
	}
	else include[l_tot] = TRUE ;
    }
    assert(l_exc == L_exc) ;
    
    FREE(exclude) ;
    return include ;
}

double *make_freqs_vector(char *freqs_file, int L_tot, int L_inc, bool *include) {
    double *freqs_tot, *freqs_inc ;
    int l_inc, l_tot ;

    freqs_tot = read_matrix_double(freqs_file, L_tot, 1, "%lf") ;
    
    freqs_inc = ALLOC(L_inc, double) ;
    l_inc = 0 ;
    for(l_tot = 0 ; l_tot < L_tot ; l_tot++)
	if(include[l_tot])
	    freqs_inc[l_inc++] = freqs_tot[l_tot] ;

    assert(l_inc == L_inc) ;
    FREE(freqs_tot) ;
    return freqs_inc ;
}

int set_monomorphs_for_exclusion(bool *snp_include, double *freqs, int L_tot, int L_inc) {
    int l_tot, l_inc, count ;

    l_inc = count = 0 ;
    for(l_tot = 0 ; l_tot < L_tot ; l_tot++)
	if(snp_include[l_tot]) {
	    if(freqs[l_inc] == 0.0 || freqs[l_inc] == 1.0) {
		count++ ;
		snp_include[l_tot] = FALSE ;
	    }
	    l_inc++ ;
	}
    assert(l_inc == L_inc) ;
    return count ;
}


/*
  dsyevr returns with eigenvalues in *increasing* order, and eigenvecs (column major)
  in corresponding order. These functions write them in the reverse order.
*/

void write_evecs(double *evecs, FILE *fp, int n, int nvecs) {
    int i, v ;
    for(v = nvecs - 1 ; v >= 0 ; v--) {
	for(i = 0 ; i < n ; i++)
	    fprintf(fp, "%lf ", evecs[i + v*n]) ;
	fputc('\n', fp) ;
    }
}

void write_evals(double *evals, FILE *fp, int nvecs) {
    int v ;
    for(v = nvecs - 1 ; v >= 0 ; v--)
	fprintf(fp, "%lf ", evals[v]) ;
    fputc('\n', fp) ;
}



/*
  Input: 1 genotype per byte:   (NA,AA,Aa,aa) <--> (00,01,02,03)  these are hexadecimal pairs representing a byte
  Output: 4 genotypes per byte: (NA,AA,Aa,aa) <--> (00,01,10,11)  these are the two-bit sequences that represent a genotype
*/

void bit_encode(unsigned char *in, int *nin_, unsigned char *out) {
    int i, j, nin = *nin_, nout ;
    
    nout = ceil(nin / 4.0) ;
    
    for(i = 0 ; i < nin ; ++out) {
	*out = 0 ;
	for(j = 3 ; j > -1 && i++ < nin ; --j, ++in)
	    *out |= *in * (int) pow(4.0, j) ; /* I don't think I'll ever understand this line again. */
    }
}

/*
unsigned char *read_data(char *genotypes_file, int *snp_range, int n, int *geno_offset) {
    
    // This is not ready for use
    
    // ERROR("Not supposed to use this function.") ;

    unsigned char *x ;
    int L = (snp_range[1] - snp_range[0]) ;
    int l, li, i, start_byte, end_byte ;
    FILE *p = fopen(genotypes_file, "rb") ;
    const int genos_per_byte = 4 ;

#if BIT_ENCODING
    x = (unsigned char *) malloc(ceil(  n * L / genos_per_byte )) ;
    
    start_geno = (snp_range[0] * n) ;
    *geno_offset = start_geno % genos_per_byte ; // Will have to skip first 'geno_offset' genotypes
    end_geno = (snp_range[1] * n) ;
    
    start_byte = start_geno / genos_per_byte ;
    end_byte = end_geno / genos_per_byte ;
    
    fseek(p, (long) start_byte, SEEK_SET) ;
    fread(x, 1, (end_byte - start_byte), p) ;
    
#else
    x = ALLOC(n * L, unsigned char) ;

    for(l = li = 0 ; l < L ; ++l)
	for(i = 0 ; i < n ; ++i, ++li)
	    x[li] = '1' ;
#endif

    fclose(p) ;
    return x ;
}

*/
