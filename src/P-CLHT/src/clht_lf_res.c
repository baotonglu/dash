/*   
 *   File: clht_lf_res.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: lock-free cache-line hash table with resizing.
 *   clht_lf_res.c is part of ASCYLIB
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *	      	      Distributed Programming Lab (LPD), EPFL
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 * FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <stdbool.h>

#include "clht_lf_res.h"

__thread ssmem_allocator_t* clht_alloc;

#ifdef DEBUG
__thread uint32_t put_num_restarts = 0;
__thread uint32_t put_num_failed_expand = 0;
__thread uint32_t put_num_failed_on_new = 0;
#endif

#include "stdlib.h"
#include "assert.h"

    const char*
clht_type_desc()
{
    return "CLHT-LF-RESIZE";
}

    inline int
is_power_of_two (unsigned int x) 
{
    return ((x != 0) && !(x & (x - 1)));
}

    static inline
int is_odd (int x)
{
    return x & 1;
}

/** Jenkins' hash function for 64-bit integers. */
    inline uint64_t
__ac_Jenkins_hash_64(uint64_t key)
{
    key += ~(key << 32);
    key ^= (key >> 22);
    key += ~(key << 13);
    key ^= (key >> 8);
    key += (key << 3);
    key ^= (key >> 15);
    key += ~(key << 27);
    key ^= (key >> 31);
    return key;
}

static uint64_t clflush_count = 0;
static uint64_t mfence_count = 0;
static uint64_t write_latency = 0;
static uint64_t CPU_FREQ_MHZ = 2100;

static inline void cpu_pause()
{
    __asm__ volatile ("pause" ::: "memory");
}

static inline unsigned long read_tsc(void)
{
    unsigned long var;
    unsigned int hi, lo;

    asm volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    var = ((unsigned long long int) hi << 32) | lo;

    return var;
}

static inline void mfence() {
#ifdef PERFCNT
    mfence_count++;
#endif
    asm volatile("mfence":::"memory");
}

static inline void clflush(char *data, int len, bool fence)
{
    volatile char *ptr = (char *)((unsigned long)data &~(CACHE_LINE_SIZE-1));
    if (fence)
        mfence();
    for(; ptr<data+len; ptr+=CACHE_LINE_SIZE){
        unsigned long etsc = read_tsc() + (unsigned long)(write_latency*CPU_FREQ_MHZ/1000);
        asm volatile(".byte 0x66; xsaveopt %0" : "+m" (*(volatile char *)(ptr)));
        while(read_tsc() < etsc) cpu_pause();
#ifdef PERFCNT
        clflush_count++;
#endif
    }
    if (fence)
        mfence();
}

/* Create a new bucket. */
    bucket_t*
clht_bucket_create() 
{
    bucket_t* bucket = NULL;
    bucket = memalign(CACHE_LINE_SIZE, sizeof(bucket_t));
    if (bucket == NULL)
    {
        return NULL;
    }

    uint32_t j;
    for (j = 0; j < KEY_BUCKT; j++)
    {
        bucket->snapshot = 0;
        bucket->key[j] = 0;
    }

    return bucket;
}

clht_hashtable_t* clht_hashtable_create(uint64_t num_buckets);

    clht_t* 
clht_create(uint64_t num_buckets)
{
    clht_t* w = (clht_t*) memalign(CACHE_LINE_SIZE, sizeof(clht_t));
    if (w == NULL)
    {
        printf("** malloc @ hatshtalbe\n");
        return NULL;
    }

    w->ht = clht_hashtable_create(num_buckets);

    w->resize_lock = 0;
    w->gc_lock = 0;
    w->status_lock = 0;
    w->version_list = NULL;
    w->version_min = 0;
    w->ht_oldest = w->ht;

    return w;
}

    clht_hashtable_t* 
clht_hashtable_create(uint64_t num_buckets) 
{
    clht_hashtable_t* hashtable = NULL;

    if (num_buckets == 0)
    {
        return NULL;
    }

    /* Allocate the table itself. */
    hashtable = (clht_hashtable_t*) memalign(CACHE_LINE_SIZE, sizeof(clht_hashtable_t));
    if (hashtable == NULL) 
    {
        printf("** malloc @ hatshtalbe\n");
        return NULL;
    }

    /* hashtable->table = calloc(num_buckets, (sizeof(bucket_t))); */
    hashtable->table = (bucket_t*) memalign(CACHE_LINE_SIZE, num_buckets * (sizeof(bucket_t)));
    if (hashtable->table == NULL) 
    {
        printf("** alloc: hashtable->table\n"); fflush(stdout);
        free(hashtable);
        return NULL;
    }

    memset((void*) hashtable->table, 0, num_buckets * (sizeof(bucket_t)));

    uint64_t i;
    for (i = 0; i < num_buckets; i++)
    {
        uint32_t j;
        for (j = 0; j < ENTRIES_PER_BUCKET; j++)
        {
            hashtable->table[i].snapshot = 0;
            hashtable->table[i].key[j] = 0;
        }
    }

    hashtable->num_buckets = num_buckets;
    hashtable->hash = num_buckets - 1;

    hashtable->table_new = NULL;
    hashtable->table_prev = NULL;

    return hashtable;
}


/* Hash a key for a particular hash table. */
    uint64_t
clht_hash(clht_hashtable_t* hashtable, clht_addr_t key) 
{
    /* uint64_t hashval; */
    /* return __ac_Jenkins_hash_64(key) & (hashtable->hash); */
    /* return hashval % hashtable->num_buckets; */
    /* return key % hashtable->num_buckets; */
    /* return key & (hashtable->num_buckets - 1); */
    return key & (hashtable->hash);
}


    static inline clht_val_t
clht_bucket_search(bucket_t* bucket, clht_addr_t key)
{
    int i;
    for (i = 0; i < KEY_BUCKT; i++)
    {
        clht_val_t val = bucket->val[i];
#ifdef __tile__
        _mm_lfence();
#endif
        if (bucket->map[i] == MAP_VALID)
        {
            if (bucket->key[i] == key)
            {
                if (likely(bucket->val[i] == val))
                {
                    return val;
                }
                else
                {
                    return 0;
                }
            }
        }
    }
    return 0;
}


/* Retrieve a key-value entry from a hash table. */
    clht_val_t
clht_get(clht_hashtable_t* hashtable, clht_addr_t key)
{
    CLHT_GC_HT_VERSION_USED(hashtable);
    size_t bin = clht_hash(hashtable, key);
    bucket_t* bucket = hashtable->table + bin;

    return clht_bucket_search(bucket, key);
}

__thread size_t num_retry_cas1 = 0, num_retry_cas2 = 0, num_retry_cas3 = 0, num_retry_cas4 = 0, num_retry_cas5 = 0;

    void
clht_print_retry_stats()
{
    printf("#cas1: %-8zu / #cas2: %-8zu / #cas3: %-8zu / #cas4: %-8zu\n",
            num_retry_cas1, num_retry_cas2, num_retry_cas3, num_retry_cas4);
}

#define DO_LF_STATS 0

#if DO_LF_STATS == 1
#  define INC(x) x++
#else
#  define INC(x) ;
#endif


/* Insert a key-value entry into a hash table. */
    int
clht_put(clht_t* h, clht_addr_t key, clht_val_t val) 
{
    int empty_retries = 0;
retry_all:
    CLHT_CHECK_RESIZE(h);
    clht_hashtable_t* hashtable = h->ht;
    size_t bin = clht_hash(hashtable, key);
    bucket_t* bucket = hashtable->table + bin;

    int empty_index = -2;
    clht_snapshot_all_t s, s1;

retry:
    s = bucket->snapshot;
#ifdef __tile__
    _mm_lfence();
#endif

    if (clht_bucket_search(bucket, key) != 0)
    {
        if (unlikely(empty_index >= 0))
        {
            bucket->map[empty_index] = MAP_INVLD;
        }
        CLHT_NO_UPDATE();
        return false;
    }

    if (likely(empty_index < 0))
    {
        empty_index = snap_get_empty_index(s);
        if (empty_index < 0)
        {
            if (empty_retries++ >= CLHT_NO_EMPTY_SLOT_TRIES)
            {
                empty_retries = 0;
                ht_status(h, 0, 2, 0);
            }
            goto retry_all;
        }
        s1 = snap_set_map(s, empty_index, MAP_INSRT);
        if (CAS_U64(&bucket->snapshot, s, s1) != s)
        {
            empty_index = -2;
            INC(num_retry_cas1);
            goto retry;
        }

        bucket->val[empty_index] = val;
#ifdef __tile__
        _mm_sfence();
#endif
        bucket->key[empty_index] = key;
#ifdef __tile__
        _mm_sfence();
#endif
        clflush((char *)&bucket->val[empty_index], sizeof(uintptr_t), true);
    }
    else
    {
        s1 = snap_set_map(s, empty_index, MAP_INSRT);
    }

    clht_snapshot_all_t s2 = snap_set_map_and_inc_version(s1, empty_index, MAP_VALID);
    if (CAS_U64(&bucket->snapshot, s1, s2) != s1)
    {
        INC(num_retry_cas2);
        goto retry;
    } else {
        clflush((char *)&bucket->snapshot, sizeof(uintptr_t), true);
    }

    CLHT_NO_UPDATE();
    return true;
}


/* Remove a key-value entry from a hash table. */
    clht_val_t
clht_remove(clht_t* h, clht_addr_t key)
{
    CLHT_CHECK_RESIZE(h);
    clht_hashtable_t* hashtable = h->ht;
    size_t bin = clht_hash(hashtable, key);
    bucket_t* bucket = hashtable->table + bin;

    clht_snapshot_t s;

    int i;
retry:
    s.snapshot = bucket->snapshot;
#ifdef __tile__
    _mm_lfence();
#endif

    for (i = 0; i < KEY_BUCKT; i++)
    {
        if (bucket->key[i] == key && s.map[i] == MAP_VALID)
        {
            clht_val_t removed = bucket->val[i];
#ifdef __tile__
            _mm_mfence();
#endif
            clht_snapshot_all_t s1 = snap_set_map(s.snapshot, i, MAP_INVLD);
            if (CAS_U64(&bucket->snapshot, s.snapshot, s1) == s.snapshot)
            {
                CLHT_NO_UPDATE();
                return removed;
            }
            else
            {
                INC(num_retry_cas3);
                goto retry;
            }
        }
    }

    CLHT_NO_UPDATE();
    return false;
}


    static uint32_t
clht_put_seq(clht_hashtable_t* hashtable, clht_addr_t key, clht_val_t val, uint64_t bin) 
{
    volatile bucket_t* bucket = hashtable->table + bin;
    uint32_t j;
    for (j = 0; j < KEY_BUCKT; j++) 
    {
        if (bucket->key[j] == 0)
        {
            bucket->val[j] = val;
            bucket->key[j] = key;
            bucket->map[j] = MAP_VALID;
            return true;
        }
    }

    printf("[CLHT] even the new ht does not have space (bucket %zu) \n", bin);
    return false;
}

    static int
bucket_cpy(volatile bucket_t* bucket, clht_hashtable_t* ht_new)
{
    uint32_t j;
    for (j = 0; j < KEY_BUCKT; j++) 
    {
        if (bucket->map[j] == MAP_VALID)
        {
            clht_addr_t key = bucket->key[j];
            uint64_t bin = clht_hash(ht_new, key);
            clht_put_seq(ht_new, key, bucket->val[j], bin);
        }
    }

    return 1;
}


/* resizing */
    int 
ht_resize_pes(clht_t* h, int is_increase, int by)
{
    /* if (!is_increase) */
    /*   { */
    /*     return 0; */
    /*   } */
    /* is_increase = 1; */
    /* if (by > 4) */
    /*   { */
    /*     by = 4; */
    /*   } */

    if (CLHT_LOCK_RESIZE(h))
    {
        /* printf("[RESPES-%02d] got resize lock\n", clht_gc_get_id()); */
        //      ht_resize_pes(h, 1, 2);
    }
    else
    {
        /* printf("[RESPES-%02d] already locked\n", clht_gc_get_id()); */
        return 0;
    }

    /* printf("[RESPES-%02d] inc: %d / by: %d\n", clht_gc_get_id(), is_increase, by); */
    //ticks s = getticks();

    clht_hashtable_t* ht_old = h->ht;

    size_t num_buckets_new;
    if (is_increase == true)
    {
        num_buckets_new = by * ht_old->num_buckets;
    }
    else
    {
        num_buckets_new = ht_old->num_buckets / 2;
    }

    clht_hashtable_t* ht_new = clht_hashtable_create(num_buckets_new);

    size_t cur_version = ht_old->version;
    ht_old->version++;

    CLHT_GC_HT_VERSION_USED(ht_old);

    size_t version_min;
    do
    {
        version_min = clht_gc_min_version_used(h);
    }
    while(cur_version >= version_min);

    ht_new->version = cur_version + 2;

    int32_t b;
    for (b = 0; b < ht_old->num_buckets; b++)
    {
        bucket_t* bu_cur = ht_old->table + b;
        bucket_cpy(bu_cur, ht_new);
    }

    ht_new->table_prev = ht_old;

    mfence();
    clflush((char *)ht_new, sizeof(clht_hashtable_t), false);
    clflush((char *)ht_new->table, num_buckets_new * sizeof(bucket_t), false);
    mfence();

    SWAP_U64((uint64_t*) h, (uint64_t) ht_new);
    clflush((char *)h, sizeof(uintptr_t), true);
    ht_old->table_new = ht_new;

    CLHT_RLS_RESIZE(h);

    //ticks e = getticks() - s;
    //printf("[RESIZE-%02d] to #bu %7zu    | took: %13llu ti = %8.6f s\n",
    //        0, ht_new->num_buckets, (unsigned long long) e, e / 2.1e9);

    return 1;
}



    size_t
clht_size(clht_hashtable_t* hashtable)
{
    uint64_t num_buckets = hashtable->num_buckets;
    bucket_t* bucket = NULL;
    size_t size = 0;

    uint64_t bin;
    for (bin = 0; bin < num_buckets; bin++)
    {
        bucket = hashtable->table + bin;
        int i;
        for (i = 0; i < KEY_BUCKT; i++)
        {
            if (bucket->key[i] != 0  && bucket->map[i] == MAP_VALID)
            {
                size++;
            }
        }
    }
    return size;
}

    size_t
ht_status(clht_t* h, int resize_increase, int emergency_increase, int just_print)
{
    if (TRYLOCK_ACQ(&h->status_lock) && !resize_increase)
    {
        return 0;
    }

    clht_hashtable_t* hashtable = h->ht;
    uint64_t num_buckets = hashtable->num_buckets;
    volatile bucket_t* bucket = NULL;
    size_t size = 0;

    uint64_t bin;
    for (bin = 0; bin < num_buckets; bin++)
    {
        bucket = hashtable->table + bin;

        uint32_t j;
        for (j = 0; j < ENTRIES_PER_BUCKET; j++)
        {
            if (bucket->key[j] > 0 && bucket->map[j] == MAP_VALID)
            {
                size++;
            }
        }
    }

    double full_ratio = 100.0 * size / ((hashtable->num_buckets) * ENTRIES_PER_BUCKET);

    if (just_print)
    {
        printf("[STATUS-%02d] #bu: %7zu / #elems: %7zu / full%%: %8.4f%% \n",
                99, hashtable->num_buckets, size, full_ratio);
    }
    else
    {
        if (full_ratio > 0 && full_ratio < CLHT_PERC_FULL_HALVE)
        {
            //printf("[STATUS-%02d] #bu: %7zu / #elems: %7zu / full%%: %8.4f%%\n",
            //        clht_gc_get_id(), hashtable->num_buckets, size, full_ratio);
            ht_resize_pes(h, 0, 33);
        }
        else if ((full_ratio > 0 && full_ratio > CLHT_PERC_FULL_DOUBLE) || emergency_increase || resize_increase)
        {
            int inc_by = (full_ratio / CLHT_OCCUP_AFTER_RES);
            int inc_by_pow2 = pow2roundup(inc_by);

            //printf("[STATUS-%02d] #bu: %7zu / #elems: %7zu / full%%: %8.4f%%\n",
            //        clht_gc_get_id(), hashtable->num_buckets, size, full_ratio);
            if (inc_by_pow2 <= 1)
            {
                inc_by_pow2 = 2;
            }
            ht_resize_pes(h, 1, inc_by_pow2);
        }
    }

    if (!just_print)
    {
        clht_gc_collect(h);
    }

    TRYLOCK_RLS(h->status_lock);
    return size;
}

    size_t
clht_size_mem(clht_hashtable_t* h) /* in bytes */
{
    if (h == NULL)
    {
        return 0;
    }

    size_t size_tot = sizeof(clht_hashtable_t**);
    size_tot += h->num_buckets * sizeof(bucket_t);
    return size_tot;
}

    size_t
clht_size_mem_garbage(clht_hashtable_t* h) /* in bytes */
{
    if (h == NULL)
    {
        return 0;
    }

    size_t size_tot = 0;
    clht_hashtable_t* cur = h->table_prev;
    while (cur != NULL)
    {
        size_tot += clht_size_mem(cur);
        cur = cur->table_prev;
    }

    return size_tot;
}

    void
clht_print(clht_hashtable_t* hashtable)
{
    uint64_t num_buckets = hashtable->num_buckets;
    bucket_t* bucket;

    printf("Number of buckets: %zu\n", num_buckets);

    uint64_t bin;
    for (bin = 0; bin < num_buckets; bin++)
    {
        bucket = hashtable->table + bin;

        printf("[[%05zu]] ", bin);

        uint32_t j;
        do
        {
            for (j = 0; j < ENTRIES_PER_BUCKET; j++)
            {
                if (bucket->key[j])
                {
                    printf("(%-5llu/%p)-> ", (long long unsigned int) bucket->key[j], (void*) bucket->val[j]);
                }
            }
            printf(" ** -> ");
        }
        while (bucket != NULL);
        printf("\n");
    }
    fflush(stdout);
}
