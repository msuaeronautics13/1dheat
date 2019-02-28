FC=gfortran
FCFLAGS=
LDFLAGS=
OBJECTS = 1dheat.o

1dheat: $(OBJECTS)
	$(FC) $(FCFLAGS) $(LDFLAGS) $(OBJECTS) -o 1dheat

%.o : %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm *.o

install:
	cp 1dheat ~/bin
