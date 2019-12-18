#-------------------------------------------------------------------------------
# TestTecplot makefile (for linux, platform specific)
# Author: Rolf Verberg
#-------------------------------------------------------------------------------

# compiler
CC = g++

# compiler flags
CFLAGS_CC =

# linker flags: i.e. mpi libraries or options like -static
LFLAGS_CC =

# source code paths

PATH_SOURCE  = ../Source
PATH_GENERAL = ../../Source

# source code files

FILES_SERIAL=\
$(PATH_GENERAL)/clString.cpp\
$(PATH_GENERAL)/clTecplot.cpp\
$(PATH_GENERAL)/clUtils.cpp\
$(PATH_SOURCE)/testtecplot.cpp

# compilation

CFLAGS = $(CFLAGS_CC) -c
LFLAGS = $(LFLAGS_CC) -lm

COBJ   = $(FILES_SERIAL:.cpp=.o)

.SUFFIXES: .o

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

all:
	$(MAKE) clobber
	$(MAKE) serial

serial: $(COBJ)
	$(CC) $(COBJ) $(LFLAGS) -o testtecplot

clobber:
	rm -f testtecplot
	rm -f $(PATH_GENERAL)/clString.o
	rm -f $(PATH_GENERAL)/clTecplot.o
	rm -f $(PATH_GENERAL)/clUtils.o
	rm -f $(PATH_SOURCE)/testtecplot.o
