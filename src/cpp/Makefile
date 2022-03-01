FLAGS =-std=c++1z -O3 -ldl -lpthread -fopenmp

main: support
	$(CXX) $(FLAGS) sqlite3.o zstddeclib.o main.cpp -o ../../bin/mzBucket

debug: FLAGS +=-pg -g
debug: support
	$(CXX) $(FLAGS) sqlite3.o zstddeclib.o main.cpp -o ../../bin/mzBucketDebug

support: sqlite3.o zstddeclib.o
