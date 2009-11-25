include make.config

BINDIR = ../bin
IDIR = .
CPPFLAGS = -I $(IDIR)
CFLAGS = -O2
LDFLAGS = -lm -lgfortran
GENO_OBJS = io.o iogeno.o dan.o 
GENO_EXECS = sstat istat flipind flipgen flipgeno gen2geno geno2gen
LA_EXECS = project snpload genocov coveigen
EXECS = lines columns columns-split $(GENO_EXECS) $(LA_EXECS)
LA_LIBS = $(LBLAS) $(LLAPACK) $(SHELLFISH_LDIR)/libcblas.a

all:		$(EXECS)
		mv $(EXECS) $(BINDIR)

$(GENO_EXECS):	$(GENO_OBJS)
$(LA_EXECS):	$(GENO_OBJS) blas-lapack-dan.o $(LA_LIBS)

## The following obscure syntax is required because io-templates.c is #included
io.o:		io.c io-templates.c io.h
		$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:		
		rm -f *.o $(EXECS)
		cd $(BINDIR) && rm -f $(EXECS)
