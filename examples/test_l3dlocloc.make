FC=gfortran 
FFLAGS= -c -O2 -fallow-argument-mismatch -std=legacy -w

default: all

all: 
	$(FC) $(FFLAGS) ../src/laprouts3d.f
	$(FC) $(FFLAGS) ../src/yrecursion.f
	$(FC) $(FFLAGS) ../src/l3dtrans.f
	$(FC) $(FFLAGS) ../src/prini.f
	$(FC) $(FFLAGS) ../src/rotviarecur3.f
	$(FC) $(FFLAGS) test_l3dlocloc.f
	$(FC) -o int2-test-l3dlocloc test_l3dlocloc.o laprouts3d.o yrecursion.o l3dtrans.o prini.o rotviarecur3.o
	./int2-test-l3dlocloc

clean: 
	rm -f *.o int2-test-l3dlocloc
