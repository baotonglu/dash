CC = gcc
CFLAGS = -O3 -Wall
LDFLAGS = -lsspfd -lm -lrt
VER_FLAGS = -D_GNU_SOURCE

ifeq ($(VERSION),DEBUG) 
CFLAGS = -O0 -ggdb -Wall -g -fno-inline
VER_FLAGS += -DDEBUG
endif

UNAME := $(shell uname -n)

all: libsspfd.a sspfd_test

default: all

sspfd_test: libsspfd.a sspfd_test.o
	$(CC) $(VER_FLAGS) -o sspfd_test sspfd_test.o $(CFLAGS) $(LDFLAGS) -L./ 

sspfd_test.o: sspfd_test.c
	$(CC) $(VER_FLAGS) -c sspfd_test.c $(CFLAGS) -L./ 

sspfd.o: sspfd.c sspfd.h
	$(CC) $(VER_FLAGS) -c sspfd.c $(CFLAGS) -L./ 	

libsspfd.a: sspfd.o sspfd.h
	@echo Archive name = libsspfd.a
	ar -r libsspfd.a sspfd.o
	rm -f *.o	
clean:
	rm -f *.o *.a sspfd_test 
