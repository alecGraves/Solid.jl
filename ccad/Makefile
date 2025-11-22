CXX := clang++
CC = clang
AR := llvm-ar
LD=lld
PREFIX ?= /usr/local

# TARGET
CFLAGS := -Wall -Wextra -Wpedantic -fPIC -flto \
          -ffunction-sections -fdata-sections \
          -fstack-protector-all -D_FORTIFY_SOURCE=3 \
          -I . -Os -mtune=native \
          -I/usr/include/opencascade

CLIBS := -lTKernel -lTKMath -lTKPrim -lTKGeomBase -lTKGeomAlgo \
	-lTKBRep -lTKTopAlgo -lTKOffset -lTKBool -lTKShHealing \
	-lTKMesh -lTKSTEP -lTKSTEPAttr -lTKSTL -lTKIGES -lTKG3d \
	-lTKBO -lstdc++ -lm

CEXE :=  -Wl,-gc-sections,-s $(CLIBS) -fuse-ld=$(LD)
# -fuse-ld=lld

LDFLAGS := -fuse-ld=$(LD)
# -static-libstdc++

ARCH ?= $(shell arch)


ifeq ($(ARCH),x86_64)
    CFLAGS += -march=native
endif
ifeq ($(ARCH),aarch64)
    CFLAGS += -march=armv8-a
endif

all: build/libccad.so build/libccad.a

build/libccad.so: build/ccad.o build/ccad_c.o
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

build/libccad.a: build/ccad.o build/ccad_c.o
	$(AR) rcs $@ $^

build/ccad.o: ccad.cpp
	mkdir -p build
	$(CXX) -c $(CFLAGS) $^ -o $@

build/ccad_c.o: ccad_c.cpp
	mkdir -p build
	$(CXX) -c $(CFLAGS) $^ -o $@

build/test_box: test/test_box.c build/libccad.a
	$(CC) -o $@ $^ -std=c11 $(CFLAGS) $(CEXE) 

build/test_loft: test/test_loft.c build/libccad.a
	$(CC) -o $@ $^ -std=c11 $(CFLAGS) $(CEXE) 

build/test_basic: test/test_basic.c build/libccad.a
	$(CC) -o $@ $^ -std=c11 $(CFLAGS) $(CEXE) 


compile_commands.json: Makefile ccad.cpp ccad_c.cpp test/test_box.c test/test_loft.c test/test_basic.c
	@echo '[\n' \
	'{"directory":"$(PWD)","command":"$(CXX) -c $(CFLAGS) ccad.cpp -o build/ccad.o","file":"ccad.cpp"},\n' \
	'{"directory":"$(PWD)","command":"$(CXX) -c $(CFLAGS) ccad_c.cpp -o build/ccad_c.o","file":"ccad_c.cpp"},\n' \
	'{"directory":"$(PWD)","command":"$(CXX) -shared -o build/libccad.so build/ccad.o build/ccad_c.o $(LDFLAGS)","file":"build/libccad.so"},\n' \
	'{"directory":"$(PWD)","command":"$(AR) rcs build/libccad.a build/ccad.o build/ccad_c.o","file":"build/libccad.a"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) -std=c11 $(CFLAGS) $(CEXE) -o build/test_box test/test_box.c build/libccad.a","file":"test/test_box.c"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) -std=c11 $(CFLAGS) $(CEXE) -o build/test_loft test/test_loft.c build/libccad.a","file":"test/test_loft.c"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) -std=c11 $(CFLAGS) $(CEXE) -o build/test_basic test/test_basic.c build/libccad.a","file":"test/test_basic.c"}\n]' \
	> $@

install:
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	mkdir -p $(PREFIX)/lib/pkgconfig
	install -m 0644 build/libccad.so $(PREFIX)/lib/
	install -m 0644 build/libccad.a $(PREFIX)/lib/
	install -m 0644 ccad.h $(PREFIX)/include/
	install -m 0644 ccad.pc $(PREFIX)/lib/pkgconfig
	sudo ldconfig

uninstall:
	rm -f $(PREFIX)/lib/libccad.so
	rm -f $(PREFIX)/lib/libccad.a
	rm -f $(PREFIX)/include/ccad.h
	rm -f $(PREFIX)/lib/pkgconfig/ccad.pc
	sudo ldconfig

format:
	clang-format -style=file -i $$(find . -type f \( -name '*.c' -o -name '*.h' -o -name '*.cpp' -o -name '*.hpp' \))

clean:
	rm -rf build
	rm -f loft.stl
	rm -f box.iges box.step box.stl box.obj

check: build/test_basic build/test_box build/test_loft
	./build/test_box
	./build/test_loft
	./build/test_basic

.PHONY: all clean check install uninstall format