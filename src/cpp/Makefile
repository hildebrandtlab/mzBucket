FLAGS =-std=c++1z -O3 -ldl -lpthread -fopenmp

main: support
	$(CXX) -std=c++1z sqlite3.o zstddeclib.o main.cpp -o ../../bin/mzBucket -O3 -ldl -lpthread -fopenmp

debug: FLAGS +=-pg -g
debug: support
	$(CXX) -std=c++1z sqlite3.o -pg -g zstddeclib.o main.cpp -o ../../bin/mzBucketDebug -O3 -ldl -lpthread -fopenmp

support: sqlite3.o zstddeclib.o
