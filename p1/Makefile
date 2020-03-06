C = icc
CFLAGS = -Wall -Werror -std=c99

VECFLAGS=-march=core-avx2 -O2 -qopt-report=5 -qopt-report-phase=vec -D NOALIAS -D ALIGNED
NOVECFLAGS=-march=core-avx2 -O2 -no-vec -D NOALIAS -D ALIGNED

all: mmvec mmnovec

clean:
	rm mmvec mmnovec mm.optrpt

mmvec: mm.c
	$(C) -o mmvec $(CFLAGS) $(VECFLAGS) mm.c

mmnovec: mm.c
	$(C) -o mmnovec $(CFLAGS) $(NOVECFLAGS) mm.c


