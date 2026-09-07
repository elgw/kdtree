CC=gcc

CFLAGS=-Wall -Wextra -pedantic  -Iinclude/
LDFLAGS=-lm -lpthread

GSL?=0

ifeq ($(GSL),1)
CFLAGS+=`pkg-config gsl --cflags` -DGSL
LDFLAGS+=`pkg-config gsl --libs`
endif

SAN?=0
DEBUG?=0

ifeq ($(SAN),1)
CFLAGS+=-fsanitize=address
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

SRC=src/kdtree.c src/pqheap.c src/quickselect.c
OBJ=kdtree.o pqheap.o quickselect.o

kdtree_ut: kdtree.o test/kdtree_ut.c makefile
	$(CC) $(CFLAGS)  test/kdtree_ut.c $(OBJ) $(LDFLAGS) -o kdtree_ut

kdtree.o: src/kdtree.c
	$(CC) -c $(CFLAGS) src/kdtree.c --std=c99

pqheap.o: src/pqheap.c
	$(CC) -c $(CFLAGS) src/pqheap.o --std=c99

quickselect.o: src/quickselect.c
	$(CC) -c $(CFLAGS) src/quickselect.o --std=c99

libkdtree.a: $(SRCFILES) makefile
	$(CC) -c $(CFLAGS) $(SRC) $(LDFLAGS)
	ar rcs libkdtree.a *.o


libkdtree.so: $(SRCFILES) makefile
	$(CC) $(CFLAGS) -fPIC -shared $(SRC) $(LDFLAGS) -o libkdtree.so

install: libkdtree.so include/kdtree.h
	# TODO COPY to DESIRED PATHS
