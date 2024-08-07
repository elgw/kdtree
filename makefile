CC=gcc

CFLAGS=-Wall -Wextra -pedantic -std=gnu11 `pkg-config gsl --cflags`
LDFLAGS=-lm -lpthread `pkg-config gsl --libs`

DEBUG?=0

ifeq ($(DEBUG),1)
CFLAGS+=-g3
else
CFLAGS+=-O3
LDFLAGS+=-flto
endif

FANALYZER?=0
ifeq ($(FANALYZER),1)
CFLAGS+=-fanalyzer
endif

kdtree_ut: src/kdtree.c makefile src/kdtree_ut.c src/pqheap.c
	$(CC) $(CFLAGS) src/kdtree.c src/kdtree_ut.c src/pqheap.c $(LDFLAGS) -o kdtree_ut


kdtree: src/kdtree.c include/kdtree.h makefile
	$(CC) $(CFLAGS) -c src/kdtree.c $(LDFLAGS)


median5_ut: src/median5.c src/median5_ut.c
	$(CC) $(CFLAGS) src/median5.c src/median5_ut.c $(LDFLAGS) -o median5

amedian: src/amedian.c
	$(CC) $(CFLAGS) src/amedian.c $(LDFLAGS) -o amedian
