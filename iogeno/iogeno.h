void read_chiamo(FILE *f, int n, int L_tot, bool *snp_include, double *dest) ;
void read_submatrix_genotypes_double(FILE *f, int majordim, bool *majorinc, int minordim, bool *minorinc, double *dest) ;
void read_submatrix_genotypes_double_transpose(FILE *f,
					       int majordim, bool *majorinc, int majorincn,
					       int minordim, bool *minorinc, int minorincn,
					       double *dest) ;
void centre_and_zero_missing_and_scale(double *g, int L, int n, double *p) ;
void read_individual_alone(char *fname, int L_tot, int L_inc, bool *include, double *dest) ;
unsigned char *read_genotypes_snp_range_all_indivs(char *file, int n, int snp_start, int snp_stop) ;
double *read_genotypes(char **fnames,
		       int L_tot, int L_inc, bool *snp_include,
		       int n_tot, int n_inc, bool *indiv_include, bool indiv_major_order) ;
int *read_individual(int i, FILE *f, int n, int L_tot, int L_inc, bool *include) ;
bool *make_include_vector(char *exclude_file, int L_tot, int L_inc) ;
double *make_freqs_vector(char *freqs_file, int L_tot, int L_inc, bool *include) ;
int set_monomorphs_for_exclusion(bool *snp_include, double *freqs, int L_tot, int L_inc) ;
void write_evecs(double *evecs, FILE *fp, int n, int nvecs) ;
void write_evals(double *evals, FILE *fp, int nvecs) ;
void bit_encode(unsigned char *in, int *nin_, unsigned char *out) ;
