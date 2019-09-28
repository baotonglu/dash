#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <thread>
#include <atomic>
#include "tbb/tbb.h"

struct timeval tv1, tv2;

using namespace std;

#include "clht_lb_res.h"
#include "ssmem.h"

void run(char **argv) {
    std::cout << "Simple Example of P-CLHT" << std::endl;

    uint64_t n = std::atoll(argv[1]);
    uint64_t *keys = new uint64_t[n];

    // Generate keys
    for (uint64_t i = 0; i < n; i++) {
        keys[i] = i + 1;
    }

    uint64_t num_thread = atoi(argv[2]);

    printf("operation,n,ops/s\n");
    typedef struct thread_data {
        uint32_t id;
        clht_t *ht;
    } thread_data_t;

    clht_t *hashtable = clht_create(512);

    thread_data_t *tds = (thread_data_t *) malloc(num_thread * sizeof(thread_data_t));

    std::atomic<int> next_thread_id;

    {
        // Load
        auto starttime = std::chrono::system_clock::now();
        next_thread_id.store(0);
        auto func = [&]() {
            int thread_id = next_thread_id.fetch_add(1);
            tds[thread_id].id = thread_id;
            tds[thread_id].ht = hashtable;

            uint64_t start_key = n / num_thread * (uint64_t)thread_id;
            uint64_t end_key = start_key + n / num_thread;

            clht_gc_thread_init(tds[thread_id].ht, tds[thread_id].id);
            ssmem_allocator_t *alloc = (ssmem_allocator_t *) malloc(sizeof(ssmem_allocator_t));
            ssmem_alloc_init_fs_size(alloc, SSMEM_DEFAULT_MEM_SIZE, SSMEM_GC_FREE_SET_SIZE, tds[thread_id].id);

            for (uint64_t i = start_key; i < end_key; i++) {
                clht_put(tds[thread_id].ht, keys[i], keys[i]);
            }
        };

        std::vector<std::thread> thread_group;

        for (int i = 0; i < num_thread; i++)
            thread_group.push_back(std::thread{func});

        for (int i = 0; i < num_thread; i++)
            thread_group[i].join();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::system_clock::now() - starttime);
        printf("Throughput: load, %f ,ops/us\n", (n * 1.0) / duration.count());
    }

    {
        // Run
        auto starttime = std::chrono::system_clock::now();
        next_thread_id.store(0);
        auto func = [&]() {
            int thread_id = next_thread_id.fetch_add(1);
            tds[thread_id].id = thread_id;
            tds[thread_id].ht = hashtable;

            uint64_t start_key = n / num_thread * (uint64_t)thread_id;
            uint64_t end_key = start_key + n / num_thread;

            for (uint64_t i = start_key; i < end_key; i++) {
                    uintptr_t val = clht_get(tds[thread_id].ht->ht, keys[i]);
                    if (val != keys[i]) {
                        std::cout << "[CLHT] wrong key read: " << val << "expected: " << keys[i] << std::endl;
                        exit(1);
                    }
            }
        };

        std::vector<std::thread> thread_group;

        for (int i = 0; i < num_thread; i++)
            thread_group.push_back(std::thread{func});

        for (int i = 0; i < num_thread; i++)
            thread_group[i].join();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::system_clock::now() - starttime);
        printf("Throughput: run, %f ,ops/us\n", (n * 1.0) / duration.count());
        double dur = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
             (double)(tv2.tv_sec - tv1.tv_sec);
    }

    delete[] keys;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("usage: %s [n] [nthreads]\nn: number of keys (integer)\nnthreads: number of threads (integer)\n", argv[0]);
        return 1;
    }

    run(argv);
    return 0;
}