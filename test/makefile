# To run, first
# git clone https://github.com/jtsiomb/kdtree.git
# and build the library

CC=gcc

CFLAGS=-Wall -Wextra -pedantic  -Iinclude/
LDFLAGS=-lm -lpthread

SAN?=0
DEBUG?=0

ifeq ($(SAN),1)
CFLAGS+=-fsanitize=adress
DEBUG=1
endif

ifeq ($(DEBUG),1)
CFLAGS+=-g3 -Og
else
CFLAGS+=-O3 -DNDEBUG
CFLAGS+=-fno-math-errno # Fine since we don't catch them in any case
CFLAGS+=-ffinite-math-only # Yes, points should be in the domain
CFLAGS+=-fno-signed-zeros # Don't care about that
CFLAGS+=-fno-trapping-math # neither this
LDFLAGS+=-flto
endif

FANALYZER?=0
ifeq ($(FANALYZER),1)
CFLAGS+=-fanalyzer
endif

LDFLAGS += -Lkdtree/ -lkdtree -Wl,-rpath='$$ORIGIN'/kdtree

jtsiomb_kdtree_test: jtsiomb_kdtree_test.c
	$(CC) $(CFLAGS) jtsiomb_kdtree_test.c $(LDFLAGS) -o jtsiomb_kdtree_test
