CC = g++
CFLAGS = -std=c++17 -march=native -Wall -O3
INCLUDES = $(shell find -name '*.hpp')
XXHASHPATH = bagminhash/xxhash/libxxhash.a

run: main
	./main

main: main.o    
	$(CC) -o main main.o $(XXHASHPATH) $(CFLAGS)

main.o: main.cpp $(INCLUDES)
	$(CC) -o main.o -c main.cpp $(CFLAGS)

test: tests
	./tests

tests: tests.o tests-main.o
	$(CC) -o tests tests.o tests-main.o $(XXHASHPATH) $(CFLAGS)

tests.o: tests.cpp tests-main.o $(INCLUDES)
	$(CC) -c tests.cpp -o tests.o $(CFLAGS)

tests-main.o: tests-main.cpp catch.hpp
	$(CC) -c tests-main.cpp -o tests-main.o $(CFLAGS)

clean:
	rm -rf *.o
