COMMON = -Wall -O3 -Wno-unused-result -faligned-new
CLANGFLAGS = -std=c++11
GPPFLAGS = -march=native -fvect-cost-model=unlimited

CC = g++
CFLAGS = $(COMMON) $(GPPFLAGS)

all: gridCHP++

gcc: CC = g++
gcc: CFLAGS = $(COMMON) $(GPPFLAGS)
gcc: all

clang: CC = clang++
clang: CFLAGS = $(COMMON) $(CLANGFLAGS)
clang: all

gridCHP++: gridCHP++.o clifford.o clustergen.o
	$(CC) $(CFLAGS) -o gridCHP++ gridCHP++.o clifford.o clustergen.o

gridCHP++.o: gridCHP++.cpp clifford.hpp clustergen.hpp
	$(CC) $(CFLAGS) -c gridCHP++.cpp

clustergen.o: clustergen.cpp clifford.hpp
	$(CC) $(CFLAGS) -c clustergen.cpp

clifford.o: clifford.cpp clifford.hpp
	$(CC) $(CFLAGS) -c clifford.cpp

clean:
	rm -f gridCHP++ *.o

.PHONY: gcc clang all clean
