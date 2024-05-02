CC = gcc
RM = /bin/rm -f
LD = ld


CFLAGS =  -O2 -Wall -Wunused
IFLAGS = -I. 
LFLAGS = -lm
ALL_FLAG = $(CFLAGS) $(IFLAGS) $(LFLAGS)


objects = util.o
          
					

all : anchor 

anchor : anchor.c
	cc anchor.c -o anchor $(ALL_FLAG)


clean : 
	$(RM) *.o anchor
	
