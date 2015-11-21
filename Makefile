.PHONY: print_itgl all clean

COMMON=Makefile src/common.hpp src/aces.hpp src/hartree_fock.hpp src/matrix.hpp src/utils.hpp src/mbpt.hpp src/ci.hpp

bin/calc:$(COMMON) src/main.cpp
	g++ src/main.cpp -O2 -std=c++11 -Ilib -o bin/calc

all:bin/calc

clean:
	rm -rf bin
	mkdir bin
