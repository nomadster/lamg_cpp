CC= g++
CFLAGS = -ansi -O3 -Wall -pedantic -g
WD = /l/disc1/home/lensi/lamg/
LIBDIR = $(WD)lib
LIBS = -L$(LIBDIR) -lboost_system -lm
INCLUDE = $(WD)include


.PHONY: clean all

all: 
	-@echo -e "Making test\n"
	-make eseguimi

clean:
	-rm -f *.o *~

LAMGLSSolver.o: LAMGLSSolver.cpp LAMGLSSolver.h Settings.h MDefs.h  OPTtypes.h MCFLSSolver.h Levels.h MtxOps.h
	$(CC) $(CFLAGS) -I$(INCLUDE) $(LIBS) -c $<

Levels.o: Levels.cpp MDefs.h OPTtypes.h Settings.h Levels.h MtxOps.h
	$(CC) $(CFLAGS) -I$(INCLUDE) $(LIBS) -c $<

IPClass.o: IPClass.C IPClass.h OPTtypes.h OPTvect.h
	$(CC) $(CFLAGS) -I$(INCLUDE) $(LIBS) -c $<

MCFInPnt.o: MCFInPnt.C MCFInPnt.h MCFClass.h OPTtypes.h IPClass.h MCFLSSolver.h OPTvect.h
	$(CC) $(CFLAGS) -I$(INCLUDE) $(LIBS) -c $<

PCGLSSlv.o: PCGLSSlv.C PCGLSSlv.h OPTtypes.h MCFLSSolver.h OPTvect.h
	$(CC) $(CFLAGS) -I$(INCLUDE) $(LIBS) -c $<

Main.o: Main.C MCFInPnt.h MCFClass.h OPTtypes.h IPClass.h MCFLSSolver.h PCGLSSlv.h OPTvect.h
	$(CC) $(CFLAGS) -I$(INCLUDE) $(LIBS) -c $<

OBJ = LAMGLSSolver.o Levels.o IPClass.o MCFInPnt.o PCGLSSlv.o Main.o

eseguimi: $(OBJ)
	$(CC) -I$(INCLUDE) $(CFLAGS) $(OBJ) $(LIBS) -o CIAO
