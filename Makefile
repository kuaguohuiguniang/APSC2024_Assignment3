CXX       = mpic++
CXXFLAGS ?= 
CPPFLAGS ?= -fopenmp -I./muparser/include -I./json/include

LDFLAGS ?= -L./muparser/lib
LIBS    ?= 

EXEC = main

SOURCES = main.cpp Jacobi.cpp

.PHONY = all $(EXEC) clean distclean

all: $(EXEC)

$(EXEC): $(SOURCES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

clean:
	$(RM) $(EXEC)

distclean: clean
	$(RM) *.csv *.out *.bak *~
