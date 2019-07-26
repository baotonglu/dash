CXX := g++
CFLAGS := -std=c++17 -march=native -I./ -lrt -pthread -O3 -g 
#CFLAGS := -I./ -lrt -lpthread -O2 -g 

Test_Cuckoo: Test_Cuckoo_O
	$(CXX) $(CFLAGS) -o src/test_cuckoo src/test_cuckoo.o

Test_Cuckoo_O: src/test_cuckoo.cpp
	$(CXX) $(CFLAGS) -c src/test_cuckoo.cpp -o src/test_cuckoo.o 

clean:
	rm -rf src/*.o
