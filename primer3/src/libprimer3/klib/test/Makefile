CC=gcc
CXX=g++
CFLAGS=-g -Wall -O2 -I..
CXXFLAGS=$(CFLAGS)
PROGS=khash_keith khash_keith2 khash_test

all:$(PROGS)

clean:
		rm -fr $(PROGS) *.dSYM a.out

khash_test:khash_test.c ../khash.h
		$(CC) $(CFLAGS) -o $@ khash_test.c
