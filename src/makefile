# Nov 2010 
#	Wai-ki Yip
#	makefile for FBAT
# modified 2021 march by Julian Hecker

# compiler, compile and link flags
#
CC = g++ 

NEWMAT10DIR = . 
INCLUDELIB  = /usr/local/include/

IFLAGS	    = -I$(NEWMAT10DIR) -I$(INCLUDELIB)
LFLAGS      = -L$(NEWMAT10DIR) -L/usr/local/lib -lnewmat

CFLAGS      = -O3 $(IFLAGS) -w  

all:	    fbat .depend

clean:	
	rm -f *.o
	rm -f .depend
	rm -f fbat
	rm -f *~

dep .depend:
	$(CC) $(CFLAGS) -MM *.cpp >> .depend

.cpp.o:
	$(CC) $(CFLAGS) -c $<

infiles    = betabase.o breakpoint.o data.o fbat.o gammabas.o main.o myutil.o normbase.o pedigree.o rand.o shell.o sufficient_stat.o xfunc.o ires.o

fbat:	$(infiles)
	$(CC) $(CFLAGS) $(LFLAGS) -o fbat $(infiles) libnewmat.a 

	
