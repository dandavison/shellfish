# You may need to change some entries in the initial part of this
# file. Comment out lines using a '#' character at the start of the
# line, and uncomment them by removing any such initial comment
# characters. Read on.


# ---------------------------------------------------
# Choose 64- or 32-bit version of CBLAS library (These are bundled
# with shellfish). The default is set to 64 bit. If you are on a 32
# bit machine change 'lib64' to just 'lib'.
LDIR = lib64

# ---------------------------------------------------
# Location of BLAS library
# If you're lucky, the following will do the trick and you won't have
# to change anything.
LBLAS = -lblas
# However, on a machine running Ubuntu 9.04 server edition I had
# specify the location of the libraryu file explicitly, by using the
# following line instead.
#  LBLAS = /usr/lib/libblas.so.3gf
# If you're having problems, search for files whose name starts with
# 'libblas' in system library directories like /usr/lib, /usr/lib/64,
# /usr/local/lib etc.

# ---------------------------------------------------
# Location of LAPACK library
LLAPACK = -llapack
# On a machine running Ubuntu 9.04 server edition I had to change that to
# LLAPACK = /usr/lib/liblapack.so.3gf
# See above for tips. You're looking for a file whose name starts with
# 'liblapack'

# ---------------------------------------------------
#
# END OF USER CONFIGURATION SECTION 
#
# Users should not have to alter the rest of this file

IDIR = .
CPPFLAGS = -I $(IDIR)

CFLAGS = -O2
LDFLAGS = -lm

GENO_OBJS = io.o iogeno.o dan.o 
GENO_EXECS = sstat istat flipind flipgen flipgeno gen2geno geno2gen
LA_EXECS = project snpload genocov coveigen
EXECS = lines columns columns-split $(GENO_EXECS) $(LA_EXECS)
LA_LIBS = $(LBLAS) $(LLAPACK) $(LDIR)/libcblas.a # -lg2c

all:		$(EXECS)
		cp $(EXECS) ..

$(GENO_EXECS):	$(GENO_OBJS)
$(LA_EXECS):	$(GENO_OBJS) blas-lapack-dan.o $(LA_LIBS)

## The following obscure syntax is required because io-templates.c is #included
io.o:		io.c io-templates.c io.h
		$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:		
		rm -f *.o $(EXECS)
		cd .. && rm -f $(EXECS)
