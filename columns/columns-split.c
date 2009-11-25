#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>

void ERROR (char *msg) {
    fputs(msg, stderr) ;
    fputc('\n', stderr) ;
    exit(2) ;
}

void usage() {
    ERROR("column-split -n input_numlines -f colnumfile -o output_dir < inputfile") ;
}

#define false 0
#define true 1
#define bool int

int main(int argc, char *argv[])
{
    
    size_t maxlinelength=0, maxidx=0;
    char *line=NULL, c, *outdir=NULL;
    FILE *infile=stdin, *cfile=NULL, *fp;
    bool invert=false, quiet=false, output_all=true;
    int *idx=NULL, *idxp, i, jin, jout, n, nin, nout, input_nlines=-1;
    char **output, **outputp, buff[1000];

    while((c = getopt(argc, argv, "f:n:o:qv")) != -1) {
	switch(c) {
	case 'f':
	    cfile = fopen(optarg, "r") ;
	    output_all = false;
	    break ;
	case 'o':
	    outdir = optarg;
	    break ;
	case 'n':
	    input_nlines = atoi(optarg);
	    break ;
	case 'q':
	    quiet = true;
	    break ;
	case 'v':
	    invert = true;
	    break ;
	case '?':
	    usage() ;
	}
    }
    
    if(input_nlines < 0 || outdir == NULL)
	usage();
    if(output_all) // implement output_all as invert with no requested columns
	invert = true;

    maxidx = 100;
    nout = 0;
    if(!output_all)
    {
	idx = (int *) malloc(maxidx * sizeof(int));
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
    }
   
    if(!quiet) fprintf(stderr, "Reading input into memory\n");

    for(i = 0 ; i < input_nlines ; ++i)
    {
	n = getline(&line, &maxlinelength, infile) - 1 ; // newline char
	if(n < 0 )
	    fprintf(stderr, "Reached end of input unexpectedly (after %d lines)", i);

	if (i == 0) // We're on the first line of input
	{
	    nin = n;
	    nout = invert ? nin - nout : nout;
	}
	else if(n != nin)
	{
	    fprintf(stderr, "Line 1 had length %d; line %d has length %d", nin, i+1, n);
	    ERROR("Input lines should be same length");
	}

	if(i == 0)
	{
	    output = (char **) malloc(nout * sizeof(char *));
	    for(jout = 0; jout < nout; ++jout)
		output[jout] = malloc(input_nlines * sizeof(char));
	}
	
	/* Now store the requested columns from the current line. */
	idxp = idx;
	jout = 0;
	outputp = output;
	for(jin = 0 ; jin < nin ; ++jin)
	{
	    if( output_all || (jin != *idxp) == invert )
		//output[jout++][i] = line[jin];
		(*outputp++)[i] = line[jin];
	    if(!output_all && jin == *idxp) ++idxp;
	}

	if (!quiet && i % 100 == 0)
	    fprintf(stderr, "%d\r", i);
    }
    
    if(!quiet) fprintf(stderr, "Writing output to files in %s\n", outdir);
    for (jout = 0; jout < nout; ++jout)
    {
	sprintf(buff, "%s/%05d", outdir, jout+1);
        fp = fopen(buff, "w");
	fwrite(output[jout], sizeof(char), input_nlines, fp);
	putc('\n', fp);
	fclose(fp);

	if (!quiet && jout % 100 == 0)
	    fprintf(stderr, "%d\r", jout);
    }

    return 0;
}
