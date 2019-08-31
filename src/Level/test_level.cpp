#include "src/Level/level.h"
#include "util/file_access.h"
#include <sstream>
#include <cstring>
#include <assert.h>
#include <thread>
#include <ctime>
#include <sys/time.h>

#define LOG(msg) std::cout << msg << "\n"
#define LAYOUT "_level"
const uint64_t POOLSIZE = (uint64_t)1024*1024*1024*10;

PMEMobjpool *pop;
LevelHashing *level;
struct timeval tv1, tv2;

struct range{
	int begin;
	int end;
};

struct my_root{
	PMEMoid _level;
};

void concurr_insert(struct range *_range){
	std::stringstream ss;
	int begin = _range->begin;
	int end = _range->end;
	Key_t key;
	std::string val;
	uint64_t offset;
	char array[64]; 
	Value_t value = (Value_t)array;
	for (int i = begin; i < end; ++i)
	{
		/* code */
		key = i;
		level->Insert(pop, key, value);
	}
}

void concurr_get(struct range *_range){
	Key_t key;
	int not_found = 0;
	for (int i = _range->begin; i < _range->end; ++i)
	{
		key = i;
		if (level->Get(pop, key) == NONE)
		{
			not_found++;
		}
	}
	printf("not found = %d\n", not_found);
}

void concurr_delete(struct range *_range){
	Key_t key;
	int not_found = 0;
	for (int i = _range->begin; i < _range->end; ++i)
	{
		key = i;
		if (level->Delete(pop, key) == false)
		{
			not_found++;
		}
	}
	printf("not found = %d\n", not_found);
}

int main(int argc, char const *argv[])
{
	assert(argc >= 4);
	int initCap = atoi(argv[1]);
	int insert_num = atoi(argv[2]);
	int thread_num = atoi(argv[3]);

	std::cout<<"The levels is "<<initCap<<std::endl;
	std::cout<<"The inserted number is "<<insert_num<<std::endl;
	std::cout<<"The thread number is "<<thread_num<<std::endl;
	//const char *file = "/mnt/pmem0/pmem_level.data";
	const char *file = "pmem_level.data";
	PMEMoid root;
	struct my_root *rr;

	if (file_exists(file) != 0)
		{
			//pop = pool<my_root>::create(file, LAYOUT, POOLSIZE, CREATE_MODE_RW);
			pop = pmemobj_create(file, LAYOUT, POOLSIZE, 0666);
			if (pop == NULL)
			{
				perror("pmemobj_create error");
				return 1;
			}
			root = pmemobj_root(pop, sizeof(struct my_root));
			rr = (my_root*)pmemobj_direct(root);
			//pmemobj_alloc(pop, &rr->_cceh, sizeof(CCEH), CCEH_TYPE, create_CCEH, &initCap);
			pmemobj_zalloc(pop, &rr->_level, sizeof(LevelHashing), LEVEL_TYPE);

			level = (LevelHashing*)pmemobj_direct(rr->_level);
			initialize_level(pop, level, &initCap);
			//Initialize_CCEH(pop, eh, initCap);

			std::cout<<"Successfully create a pool"<<std::endl;
		}else{
			
			pop = pmemobj_open(file, LAYOUT);
			if (pop == NULL)
			{
				perror("pmemobj_open error");
				return 1;
			}

			root = pmemobj_root(pop, sizeof(struct my_root));
			rr = (struct my_root*)pmemobj_direct(root);
			level = (LevelHashing*)pmemobj_direct(rr->_level);
			remapping(level);
			std::cout<<"Successfully open a pool"<<std::endl;
		}
		level->display_size();

int chunk_size = insert_num/thread_num;
	struct range rarray[thread_num];
	for (int i = 0; i < thread_num; ++i)
	{

		rarray[i].begin = i*chunk_size + 1;
		rarray[i].end = (i+1)*chunk_size + 1;
	}
	rarray[thread_num-1].end = insert_num + 1;

//-----------------------------------------------Concurrent Insertion Test-----------------------------------------------------------------------
	std::thread *thread_array[thread_num];
	
	LOG("Concurrent insertion begin");
	gettimeofday(&tv1, NULL);
	for (int i = 0; i < thread_num; ++i)
	{
		thread_array[i] = new std::thread(concurr_insert, &rarray[i]);
	}
	// Join operation
	for (int i = 0; i < thread_num; ++i)
	{
		// code 
		thread_array[i]->join();
		delete thread_array[i];
	}
	gettimeofday(&tv2, NULL);
	auto duration = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
		printf ("For %d threads, Insertion Total time = %f seconds, the throughput is %f options/s\n", thread_num,
	         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
	         (double) (tv2.tv_sec - tv1.tv_sec), insert_num/duration);	
	//eh->Get_Number();



//-----------------------------------------------Concurrent postive Get Test-----------------------------------------------------------------------
	//std::cout<<"There are "<<eh->GetItemNum()<<" items inserted in the hashing index!"<<std::endl;
	//rarray[thread_num-1].end = insert_num + 5;
	LOG("Concurrent positive get begin!");
	gettimeofday(&tv1, NULL);
	for (int i = 0; i < thread_num; ++i)
	{
		thread_array[i] = new std::thread(concurr_get, &rarray[i]);
	}

	for (int i = 0; i < thread_num; ++i)
	{
		thread_array[i]->join();
		delete thread_array[i];
	}
	gettimeofday(&tv2, NULL);
	LOG("Concurrent positive get done!");
	duration = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
	printf ("For %d threads, Postive Search Total time = %f seconds, the throughput is %f options/s\n", thread_num,
	         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
	         (double) (tv2.tv_sec - tv1.tv_sec), insert_num/duration);	
	//eh->Get_Number();
	
	LOG("Concurrent delete begin!");
	gettimeofday(&tv1, NULL);
	for (int i = 0; i < thread_num; ++i)
	{
		thread_array[i] = new std::thread(concurr_delete, &rarray[i]);
	}

	for (int i = 0; i < thread_num; ++i)
	{
		thread_array[i]->join();
		delete thread_array[i];
	}
	gettimeofday(&tv2, NULL);
	LOG("Concurrent delete done!");
	duration = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
	printf ("For %d threads, Delete Total time = %f seconds, the throughput is %f options/s\n", thread_num,
	         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
	         (double) (tv2.tv_sec - tv1.tv_sec), insert_num/duration);	
//----------------------------------------------Concurrent negative Get Test-----------------------------------------------------------------------
	/*
	LOG("Concurrent negative get begin!");
	for (int i = 0; i < thread_num; ++i)
	{
		rarray[i].begin = insert_num + i*chunk_size + 1;
		rarray[i].end = insert_num + (i+1)*chunk_size + 1;
	}

	gettimeofday(&tv1, NULL);
	for (int i = 0; i < thread_num; ++i)
	{
		thread_array[i] = new std::thread(concurr_get, &rarray[i]);
	}

	for (int i = 0; i < thread_num; ++i)
	{
		thread_array[i]->join();
		delete thread_array[i];
	}
	gettimeofday(&tv2, NULL);
	LOG("Concurrent negative get done!");
	duration = (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec);
	printf ("For %d threads, Negative Search Total time = %f seconds, the throughput is %f options/s\n", thread_num,
	         (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
	         (double) (tv2.tv_sec - tv1.tv_sec), insert_num/duration);	
	pmemobj_close(pop);
//-----------------------------------------------Persistence Test-----------------------------------------------------------------------
/*
	std::cout<<"Insert Testing"<<std::endl;
	pop = pmemobj_open(file, LAYOUT);
	if (pop == NULL)
	{
		perror("pmemobj_open error");
		return 1;
	}

	root = pmemobj_root(pop, sizeof(struct my_root));
	rr = (struct my_root*)pmemobj_direct(root);
	level = (LevelHashing*)pmemobj_direct(rr->_level);
	remapping(level);

	Key_t key;
	char array[64]; 
	Value_t value = (Value_t)array;

	for (int i = 1; i < insert_num+3; ++i)
	{
		key = i;
		level->Insert(pop, key, value);
	}
	pmemobj_close(pop);

	// Test if it has been correclty stored in persistent storage
	std::cout<<"Read Testing"<<std::endl;
	pop = pmemobj_open(file, LAYOUT);
	if (pop == NULL)
	{
		perror("pmemobj_open error");
		return 1;
	}

	root = pmemobj_root(pop, sizeof(struct my_root));
	rr = (struct my_root*)pmemobj_direct(root);
	level = (LevelHashing*)pmemobj_direct(rr->_level);
	remapping(level);

	for (int i = 1; i < insert_num+3; ++i)
	{
		key = i;
		if (level->Get(pop, key) == NONE)
		{
			std::cout<<"Search the key "<< i << ": ERROR!"<<std::endl;
		}
	}
	pmemobj_close(pop);
*/
	return 0;
}
