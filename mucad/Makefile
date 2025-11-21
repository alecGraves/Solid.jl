CC := clang++
AR := llvm-ar
TARGET := normal
PREFIX ?= /usr/local

CFLAGS := -Wall -Wextra -Wpedantic -Werror -fPIC -flto \
          -ffunction-sections -fdata-sections \
          -fno-stack-protector -D_FORTIFY_SOURCE=0 \
          -I . -O3 -mtune=native
CEXE := -fuse-ld=lld -Wl,-gc-sections,-icf=safe,-s
LDFLAGS := -fuse-ld=lld

ifeq ($(TARGET),embedded)
    SOURCE := mucad_embedded.c
else
    SOURCE := mucad.c
endif

ifeq ($(ARCH),x86_64)
    CFLAGS += -march=native
endif
ifeq ($(ARCH),aarch64)
    CFLAGS += -march=armv8-a
endif

all: build/libmucad.so build/libmucad.a

build/libmucad.so: build/mucad.o
	$(CC) -shared -o $@ $< $(LDFLAGS)

build/libmucad.a: build/mucad.o
	$(AR) rcs $@ $<

build/mucad.o: $(SOURCE)
	mkdir -p build
	$(CC) -c $(CFLAGS) $< -o $@

build/test: test.c build/libmucad.a
	$(CC) $(CFLAGS) $(CEXE) -o $@ test.c build/libmucad.a

build/test_all: test_all.c build/libmucad.a
	$(CC) $(CFLAGS) $(CEXE) -o $@ test_all.c build/libmucad.a

compile_commands.json: Makefile $(SOURCE)
	@echo '[\n' \
	'{"directory":"$(PWD)","command":"$(CC) -c $(CFLAGS) $(SOURCE) -o build/mucad.o","file":"$(SOURCE)"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) -shared -o build/libmucad.so build/mucad.o $(LDFLAGS)","file":"build/libmucad.so"},\n' \
	'{"directory":"$(PWD)","command":"$(AR) rcs build/libmucad.a build/mucad.o","file":"build/libmucad.a"},\n' \
	'{"directory":"$(PWD)","command":"$(CC) $(CFLAGS) $(CEXE) -o build/test test.c build/libmucad.a -fuse-ld=lld -Wl,-gc-sections,-icf=safe","file":"test.c"}\n]' \
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
	clang-format -style=file -i $$(find . -type f \( -name '*.c' -o -name '*.h' \))

clean:
	rm -rf build

check: build/test build/test_all
	./build/test
	./build/test_all

.PHONY: all clean check install uninstall format