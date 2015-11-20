.PHONY: print_itgl all clean

COMMON=Makefile src/aces.hpp src/hartree_fock.hpp src/matrix.hpp src/utils.hpp

bin/calc:$(COMMON) src/main.cpp
	g++ src/main.cpp -O -std=c++11 -Ilib -o bin/calc

all:bin/calc

clean:
	make -C print_itgl clean
	rm -rf build
	mkdir build
	rm -rf samples
	mkdir samples
	rm -rf bin
	mkdir bin
