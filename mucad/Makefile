CXX := clang++
CC=clang
AR := llvm-ar-18
PREFIX ?= /usr/local

# TARGET
CFLAGS := -Wall -Wextra -Wpedantic -fPIC -flto \
          -ffunction-sections -fdata-sections \
          -fstack-protector-all -D_FORTIFY_SOURCE=3 \
          -I . -O3 -mtune=native \
          -I/usr/include/opencascade

LIBS := -lTKernel -lTKMath -lTKPrim -lTKGeomBase -lTKGeomAlgo \
	-lTKBRep -lTKTopAlgo -lTKOffset -lTKBool -lTKShHealing \
	-lTKMesh -lTKSTEP -lTKSTEPAttr -lTKSTL -lTKIGES -lTKG3d \
	-lTKBO -lstdc++ 

CEXE := -fuse-ld=lld -Wl,-gc-sections,-icf=safe,-s $(LIBS)
LDFLAGS := -fuse-ld=lld 
# -static-libstdc++

ARCH ?= $(shell arch)


ifeq ($(ARCH),x86_64)
    CFLAGS += -march=native
endif
ifeq ($(ARCH),aarch64)
    CFLAGS += -march=armv8-a
endif

all: build/libmucad.so build/libmucad.a

build/libmucad.so: build/mucad.o build/mucad_c.o
	$(CXX) -shared -o $@ $^ $(LDFLAGS)

build/libmucad.a: build/mucad.o build/mucad_c.o
	$(AR) rcs $@ $^

build/mucad.o: mucad.cpp
	mkdir -p build
	$(CXX) -c $(CFLAGS) $^ -o $@

build/mucad_c.o: mucad_c.cpp
	mkdir -p build
	$(CXX) -c $(CFLAGS) $^ -o $@

build/test_box: test/test_box.c build/libmucad.a
	$(CC) -std=c11 $(CFLAGS) $(CEXE) -o $@ $^

build/test_loft: test/test_loft.c build/libmucad.a
	$(CC) -std=c11 $(CFLAGS) $(CEXE) -o $@ $^

build/test_basic: test/test_basic.c build/libmucad.a
	$(CC) -std=c11 $(CFLAGS) $(CEXE) -o $@ $^


compile_commands.json: Makefile mucad.cpp mucad_c.cpp test/test_box.c test/test_loft.c test/test_basic.c
	@echo '[\n' \
	'{"directory":"$(PWD)","command":"$(CXX) -c $(CFLAGS) mucad.cpp -o build/mucad.o","file":"mucad.cpp"},\n' \
	'{"directory":"$(PWD)","command":"$(CXX) -c $(CFLAGS) mucad_c.cpp -o build/mucad_c.o","file":"mucad_c.cpp"},\n' \
	'{"directory":"$(PWD)","command":"$(CXX) -shared -o build/libmucad.so build/mucad.o build/mucad_c.o $(LDFLAGS)","file":"build/libmucad.so"},\n' \
	'{"directory":"$(PWD)","command":"$(AR) rcs build/libmucad.a build/mucad.o build/mucad_c.o","file":"build/libmucad.a"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) -std=c11 $(CFLAGS) $(CEXE) -o build/test_box test/test_box.c build/libmucad.a","file":"test/test_box.c"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) -std=c11 $(CFLAGS) $(CEXE) -o build/test_loft test/test_loft.c build/libmucad.a","file":"test/test_loft.c"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) -std=c11 $(CFLAGS) $(CEXE) -o build/test_basic test/test_basic.c build/libmucad.a","file":"test/test_basic.c"}\n]' \
	> $@

install:
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include
	install -m 0644 build/libmucad.so $(PREFIX)/lib/
	install -m 0644 build/libmucad.a $(PREFIX)/lib/
	install -m 0644 mucad.h $(PREFIX)/include/
	sudo ldconfig

uninstall:
	rm -f $(PREFIX)/lib/libmucad.so
	rm -f $(PREFIX)/lib/libmucad.a
	rm -f $(PREFIX)/include/mucad.h
	sudo ldconfig

format:
	clang-format -style=file -i $$(find . -type f \( -name '*.c' -o -name '*.h' -o -name '*.cpp' -o -name '*.hpp' \))

clean:
	rm -rf build

check: build/test_basic build/test_box build/test_loft
	./build/test_box
	./build/test_loft
	./build/test_basic

.PHONY: all clean check install uninstall format