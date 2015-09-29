.PHONY: print_itgl all clean

print_itgl:
	make -C print_itgl
	mv print_itgl/xprint_itgl build

all:
	md build
	gcc --version

clean:
	make -C print_itgl clean
	rm -rf build
