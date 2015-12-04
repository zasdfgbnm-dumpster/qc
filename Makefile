.PHONY: print_itgl all clean

COMMON=Makefile src/common.hpp src/aces.hpp src/hartree_fock.hpp src/matrix.hpp src/utils.hpp src/mbpt.hpp src/ci.hpp src/tdhf.hpp

bin/calc:$(COMMON) src/main.cpp
	g++ src/main.cpp -O2 -std=c++11 -Ilib -o bin/calc

bin/debug:$(COMMON) src/main.cpp
	g++ src/main.cpp -g -O0 -std=c++11 -Ilib -o bin/debug

all:bin/calc bin/debug

clean:
	rm -rf bin
	mkdir bin
