CC = g++

ROOTLIB=$(shell root-config --libs)
ROOTINC=$(shell root-config --cflags)

FLAG =  -O3 -g -D_FILE_OFFSET_BITS=64 -std=c++11
MFILE_LIB = 
RUNFILE = r2root

INCLUDE =
LIBS = 

sourcefile=main.cpp tree.cpp GEBSort.h tree.h

all: $(RUNFILE)

$(RUNFILE): $(sourcefile)
	$(CC) $(FLAG) $(ROOTLIB) $(ROOTINC) -o $@ $(filter %.cpp, $(sourcefile))
#treeDict.cpp: GEBSort.h tree.h geb_linkdef.h
#	rootcint -f treeDict.cpp -c tree.h geb_linkdef.h

clean:
	rm -f $(RUNFILE) $(OBJFILES) 
