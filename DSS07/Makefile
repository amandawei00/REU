COMP = gfortran
FFLAGS = -O2    # This makes the code run faster! Optimizes gfortran...
RUNFLAGS = -ffixed-line-length-none -cpp -fbacktrace -fcheck=all
DEBUGFLAGS = -g -Wall -Wextra -Wsurprising -Walign-commons -ffpe-trap=zero,invalid,overflow,underflow
SOURCE = generate.f fDSS.f

executables:
	$(COMP) $(FFLAGS) -c $(SOURCE)
	$(COMP) $(LIBPATH) $(LIB) *.o -o generate.out

clean:
	-rm *.out
#	-rm *mod
	-rm *.o
#	-rm fort.*
	-rm *~
