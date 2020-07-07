
all : seg

seg : seg.c lnfac.h genwin.h genwin.o
	cc -O -o seg seg.c genwin.o -lm

hiseg : hiseg.c lnfac.h genwin.h genwin.o
	cc -O -o hiseg hiseg.c genwin.o -lm

genwin.o : genwin.c genwin.h
	cc -O -c genwin.c

clean:
	rm -f seg seg.o genwin.o

