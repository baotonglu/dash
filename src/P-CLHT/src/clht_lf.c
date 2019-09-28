/*   
 *   File: clht_lf.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: lock-free cache-line hash table with no resizing. If there is
 *    not enough space for a key/value pair in its corresponding bucket, the 
 *    operation might never complete. Thus, better use the resize version.
 *   clht_lf.c is part of ASCYLIB
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

#include "clht_lf.h"

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
    return "CLHT-LF-NO-RESIZE";
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
        return false;
    }

    if (likely(empty_index < 0))
    {
        empty_index = snap_get_empty_index(s);
        if (empty_index < 0)
        {
            goto retry;
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
    }

    return true;
}


/* Remove a key-value entry from a hash table. */
    clht_val_t
clht_remove(clht_t* h, clht_addr_t key)
{
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
                return removed;
            }
            else
            {
                INC(num_retry_cas3);
                goto retry;
            }
        }
    }
    return 0;
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


    void
clht_print(clht_hashtable_t* hashtable)
{
    uint64_t num_buckets = hashtable->num_buckets;
    bucket_t* bucket;

    printf("Number of buckets: %u\n", num_buckets);

    uint64_t bin;
    for (bin = 0; bin < num_buckets; bin++)
    {
        bucket = hashtable->table + bin;

        printf("[[%05d]] ", bin);

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
