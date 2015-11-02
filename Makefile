.PHONY: print_itgl all clean

COMMON=Makefile src/aces.hpp src/hartree_fock.hpp src/matrix.hpp src/utils.hpp

all:integrals bin/rhf

bin/rhf:$(COMMON) src/rhf.hpp src/rhf.cpp
	icpc src/rhf.cpp -O0 -qopenmp -mkl=parallel -lpthread -lm -std=c++11 -Ilib -o bin/rhf

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
