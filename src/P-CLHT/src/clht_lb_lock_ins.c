/*   
 *   File: clht_lb_lock_ins.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: lock-based cache-line hash table with no resizing. Only
 *    inserts use real locks. Removals work with a CAS.
 *   clht_lb_lock_ins.c is part of ASCYLIB
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

#include "clht_lb_lock_ins.h"

__thread ssmem_allocator_t* clht_alloc;

#ifdef DEBUG
__thread uint32_t put_num_restarts = 0;
__thread uint32_t put_num_failed_expand = 0;
__thread uint32_t put_num_failed_on_new = 0;
#endif

__thread size_t check_ht_status_steps = CLHT_STATUS_INVOK_IN;

#include "stdlib.h"
#include "assert.h"

const char*
clht_type_desc()
{
  return "CLHT-LB-LOCK-INS";
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

  bucket->lock = 0;

  uint32_t j;
  for (j = 0; j < ENTRIES_PER_BUCKET; j++)
    {
      bucket->key[j] = 0;
    }
  bucket->next = NULL;

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
  if (w->ht == NULL)
    {
      free(w);
      return NULL;
    }
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

  memset(hashtable->table, 0, num_buckets * (sizeof(bucket_t)));
    
  uint64_t i;
  for (i = 0; i < num_buckets; i++)
    {
      hashtable->table[i].lock = LOCK_FREE;
      uint32_t j;
      for (j = 0; j < ENTRIES_PER_BUCKET; j++)
	{
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


/* Retrieve a key-value entry from a hash table. */
clht_val_t
clht_get(clht_hashtable_t* hashtable, clht_addr_t key)
{
  size_t bin = clht_hash(hashtable, key);
  volatile bucket_t* bucket = hashtable->table + bin;

  uint32_t j;
  do 
    {
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) 
	{
	  clht_val_t val = bucket->val[j];
#ifdef __tile__
	  _mm_lfence();
#endif
	  if (bucket->key[j] == key) 
	    {
	      if (bucket->val[j] == val)
		{
		  return val;
		}
	      else
		{
		  return 0;
		}
	    }
	}

      bucket = bucket->next;
    } while (bucket != NULL);
  return 0;
}

static inline int
bucket_exists(volatile bucket_t* bucket, clht_addr_t key)
{
  uint32_t j;
  do 
    {
      for (j = 0; j < ENTRIES_PER_BUCKET; j++)
      	{
      	  if (bucket->key[j] == key)
      	    {
      	      return true;
      	    }
      	}
      bucket = bucket->next;
    } while (bucket != NULL);
  return false;
}

/* Insert a key-value entry into a hash table. */
int
clht_put(clht_t* h, clht_addr_t key, clht_val_t val) 
{
  clht_hashtable_t* hashtable = h->ht;
  size_t bin = clht_hash(hashtable, key);
  volatile bucket_t* bucket = hashtable->table + bin;

#if CLHT_READ_ONLY_FAIL == 1
  if (bucket_exists(bucket, key))
    {
      return false;
    }
#endif

  clht_lock_t* lock = &bucket->lock;
  LOCK_ACQ(lock, hashtable);

  volatile clht_addr_t* empty = NULL;
  clht_val_t* empty_v = NULL;

  uint32_t j;
  do 
    {
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) 
	{
	  if (bucket->key[j] == key) 
	    {
	      LOCK_RLS(lock);
	      return false;
	    }
	  else if (empty == NULL && bucket->key[j] == 0)
	    {
	      empty = &bucket->key[j];
	      empty_v = &bucket->val[j];
	    }
	}
        
      if (bucket->next == NULL)
	{
	  if (empty == NULL)
	    {
	      DPP(put_num_failed_expand);
	      bucket->next = clht_bucket_create(hashtable);
	      bucket->next->val[0] = val;
#ifdef __tile__
	      /* keep the writes in order */
	      _mm_sfence();
#endif
	      bucket->next->key[0] = key;
	    }
	  else 
	    {
	      *empty_v = val;
#ifdef __tile__
	      /* keep the writes in order */
	      _mm_sfence();
#endif
	      *empty = key;
	    }

	  LOCK_RLS(lock);
	  return true;
	}
      bucket = bucket->next;
    }
  while (true);
}


/* Remove a key-value entry from a hash table. */
clht_val_t
clht_remove(clht_t* h, clht_addr_t key)
{
  clht_hashtable_t* hashtable = h->ht;
  size_t bin = clht_hash(hashtable, key);
  volatile bucket_t* bucket = hashtable->table + bin;

  uint32_t j;
  do 
    {
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) 
	{
	  if (bucket->key[j] == key) 
	    {
	      if (CAS_U64(&bucket->key[j], key, KEY_BLCK) == key)
		{
		  clht_val_t val = bucket->val[j];
#ifdef __tile__
		  _mm_lfence();
#endif
		  bucket->key[j] = 0;
		  return val;
		}
	      else
		{
		  return false;
		}
	    }
	}
      bucket = bucket->next;
    } 
  while (bucket != NULL);
  return false;
}

static uint32_t
clht_put_seq(clht_hashtable_t* hashtable, clht_addr_t key, clht_val_t val, uint64_t bin) 
{
  volatile bucket_t* bucket = hashtable->table + bin;
  uint32_t j;

  do 
    {
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) 
	{
	  if (bucket->key[j] == 0)
	    {
	      bucket->val[j] = val;
	      bucket->key[j] = key;
	      return true;
	    }
	}
        
      if (bucket->next == NULL)
	{
	  DPP(put_num_failed_expand);
	  bucket->next = clht_bucket_create(hashtable);
	  bucket->next->val[0] = val;
	  bucket->next->key[0] = key;
	  return true;
	}

      bucket = bucket->next;
    } 
  while (true);
}


static int
bucket_cpy(volatile bucket_t* bucket, clht_hashtable_t* ht_new)
{
  if (!LOCK_ACQ_RES(&bucket->lock))
    {
      return 0;
    }
  uint32_t j;
  do 
    {
      for (j = 0; j < ENTRIES_PER_BUCKET; j++) 
	{
	  clht_addr_t key = bucket->key[j];
	  if (key != 0) 
	    {
	      uint64_t bin = clht_hash(ht_new, key);
	      clht_put_seq(ht_new, key, bucket->val[j], bin);
	    }
	}
      bucket = bucket->next;
    } 
  while (bucket != NULL);

  return 1;
}


void
ht_resize_help(clht_hashtable_t* h)
{
  if ((int32_t) FAD_U32((volatile uint32_t*) &h->is_helper) <= 0)
    {
      return;
    }

  int32_t b;
  /* hash = num_buckets - 1 */
  for (b = h->hash; b >= 0; b--)
    {
      bucket_t* bu_cur = h->table + b;
      if (!bucket_cpy(bu_cur, h->table_tmp))
	{	    /* reached a point where the resizer is handling */
	  /* printf("[GC-%02d] helped  #buckets: %10zu = %5.1f%%\n",  */
	  /* 	 clht_gc_get_id(), h->num_buckets - b, 100.0 * (h->num_buckets - b) / h->num_buckets); */
	  break;
	}
    }

  h->helper_done = 1;
}

size_t
clht_size(clht_hashtable_t* hashtable)
{
  uint64_t num_buckets = hashtable->num_buckets;
  volatile bucket_t* bucket = NULL;
  size_t size = 0;

  uint64_t bin;
  for (bin = 0; bin < num_buckets; bin++)
    {
      bucket = hashtable->table + bin;
       
      uint32_t j;
      do
	{
	  for (j = 0; j < ENTRIES_PER_BUCKET; j++)
	    {
	      if (bucket->key[j] > 0)
		{
		  size++;
		}
	    }

	  bucket = bucket->next;
	}
      while (bucket != NULL);
    }
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
  size_tot += (h->num_buckets + h->num_expands) * sizeof(bucket_t);
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
  volatile bucket_t* bucket;

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

	  bucket = bucket->next;
	  printf(" ** -> ");
	}
      while (bucket != NULL);
      printf("\n");
    }
  fflush(stdout);
}


