/*   
 *   File: prand.h
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: 
 *   prand.h is part of ASCYLIB
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

#ifndef _PRAND_H_
#define _PRAND_H_

#include <assert.h>
#include "utils.h"

#define PRAND_LEN 8192
/* #define PRAND_LEN 65536 */
/* #define PRAND_LEN 2097152 */

typedef uint32_t prand_gen_t;

static inline prand_gen_t*
prand_new()
{
  prand_gen_t* g = (prand_gen_t*) malloc(sizeof(prand_gen_t) * PRAND_LEN);
  assert(g != NULL);

  unsigned long* s = seed_rand();

  int i;
  for (i = 0; i < PRAND_LEN; i++)
    {
      g[i] = (uint32_t) xorshf96(s, s+1, s+2);
    }

  free(s);
  return g;
}

static inline prand_gen_t*
prand_new_range(size_t min, size_t max)
{
  prand_gen_t* g = (prand_gen_t*) malloc(sizeof(prand_gen_t) * PRAND_LEN);
  assert(g != NULL);

  unsigned long* s = seed_rand();

  int i;
  for (i = 0; i < PRAND_LEN; i++)
    {
      g[i] = (uint32_t) (xorshf96(s, s+1, s+2) % max) + min;
    }

  free(s);
  return g;
}

static inline prand_gen_t*
prand_new_range_len(size_t len, size_t min, size_t max)
{
  prand_gen_t* g = (prand_gen_t*) malloc(sizeof(prand_gen_t) * (len + 1));
  assert(g != NULL);

  unsigned long* s = seed_rand();

  int i;
  for (i = 0; i <= len; i++)
    {
      g[i] = (uint32_t) ((xorshf96(s, s+1, s+2) % max) + min);
      assert(g[i] != 0);
    }

  free(s);
  return g;
}


static inline prand_gen_t
prand_nxt(const prand_gen_t* g, int* idx)
{
  return g[*idx++ & (PRAND_LEN - 1)];
}

#define PRAND_GET_NXT(g, idx)			\
  g[idx++ & (PRAND_LEN - 1)]

#define PRAND_GET_NXT_L(g, idx, len)		\
  g[idx++ & (len)]


#define PRAND_FOR(g, i, key)			\
  for (i = 0; i < PRAND_LEN; key = g[i], i++)



  


#endif
