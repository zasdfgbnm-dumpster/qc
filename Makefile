.PHONY: print_itgl all clean

COMMON=Makefile src/aces.hpp src/hartree_fock.hpp src/matrix.hpp src/utils.hpp

bin/calc:$(COMMON) src/main.cpp
	#icpc src/main.cpp -O -qopenmp -mkl=parallel -lpthread -lm -std=c++11 -Ilib -o bin/calc
	module load gcc/5.2.0;\
	g++ src/main.cpp -O -std=c++11 -Ilib -o bin/calc

all:integrals bin/calc

integrals:print_itgl
	rm -rf build
	mkdir build
	cp -r task_specific/* build
	find build -mindepth 1 -maxdepth 1 -type d -exec cp bin/xprint_itgl {} \;
	for i in build/*; do \
		cp aces/* $$i;\
		cd $$i;\
		./xjoda > out;\
		./xvmol >> out;\
		./xvmol2ja >> out;\
		./xvscf >> out;\
		./xprint_itgl > $$(basename $$i.txt);\
		mv $$(basename $$i.txt) ..;\
		mv out ../$$(basename $$i.out);\
		cd ..;\
		rm -rf $$(basename $$i);\
	done
	cp build/* samples
	
print_itgl:
	module load intel/2016.0.109;\
	make -C print_itgl
	mv print_itgl/xprint_itgl bin/

clean:
	make -C print_itgl clean
	rm -rf build
	mkdir build
	rm -rf samples
	mkdir samples
	rm -rf bin
	mkdir bin
