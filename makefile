CC = g++
CFLAGS = -Wall -O2
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz -lm -lpthread #-lgmpxx -lgmp
DEBUG = -ggdb -g
SRC = src

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

pangrowth: $(SRC)/pangrowth.cpp $(SRC)/size_pangenome.h $(SRC)/khashl.h $(SRC)/ketopt.h $(SRC)/kseq.h $(SRC)/kthread.h $(SRC)/yak-hist.h
	$(CC) $(CFLAGS) -o $@ $(SRC)/pangrowth.cpp $(SRC)/kthread.c $(SRC)/size_pangenome.h $(SRC)/yak-hist.h $(LIBS)
