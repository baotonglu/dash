CXX := g++
CFLAGS := -std=c++17 -march=native -I./ -lrt -pthread -lpmemobj -O3 -g 
#CFLAGS := -I./ -lrt -lpthread -O2 -g 

Test_Cuckoo: Test_Cuckoo_O
	$(CXX) $(CFLAGS) -o src/test_cuckoo src/test_cuckoo.o

Test_Cuckoo_O: src/test_cuckoo.cpp
	$(CXX) $(CFLAGS) -c src/test_cuckoo.cpp -o src/test_cuckoo.o 

Test: Test_CCEH
	$(CXX) $(CFLAGS) -g -o src/CCEH/test_cceh src/CCEH/test_cceh.o 

Test_CCEH: src/CCEH/test_cceh.cpp src/CCEH/CCEH_base.h util/hash.h
	$(CXX) $(CFLAGS) -g -c src/CCEH/test_cceh.cpp -o src/CCEH/test_cceh.o -DINPLACE

clean:
	rm -rf src/*.o
