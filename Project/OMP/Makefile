CC=gcc
CFLAGS=-fopenmp -O3
LDFLAGS=-lm
TARGET=ballAlg.o ballQuery.o

all: $(TARGET)

ballAlg.o:
	$(CC) $(CFLAGS) $^ -o ballAlg.o $(LDFLAGS)

ballQuery.o:
	$(CC) $(CFLAGS) $^ -o ballQuery.o $(LDFLAGS)

ballAlg.o: src/ballAlg.c src/ballAlg.h src/pointArith.c src/pointArith.h lib/msort.c lib/msort.h lib/genPoints.c lib/genPoints.h
ballQuery.o: lib/ballQuery.c

clean:
	@echo Cleaning...
	rm -f $(TARGET)
