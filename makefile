IDIR = "."
CPPFLAGS = -I $(IDIR)
CFLAGS = -O2 #-Wall
LDFLAGS = -lm -llapack # -static /usr/lib64/libm.a

GENO_OBJS = io.o iogeno.o dan.o 
GENO_EXECS = sstat istat flipind flipgen flipgeno gen2geno geno2gen
LA_EXECS = project snpload genocov # coveigen
EXECS = lines columns columns-split $(GENO_EXECS) $(LA_EXECS)

all:		$(EXECS)
		cp $(EXECS) ..

$(GENO_EXECS):	$(GENO_OBJS)
$(LA_EXECS):	$(GENO_OBJS) blas-lapack-dan.o

## The following obscure syntax is required because io-templates.c is #included
io.o:		io.c io-templates.c io.h
		$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:		
		rm -f *.o $(EXECS)
		cd .. && rm -f $(EXECS)
