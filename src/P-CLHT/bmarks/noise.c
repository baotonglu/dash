/*   
 *   File: noise.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: 
 *   noise.c is part of ASCYLIB
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

#include "utils.h"
#include <assert.h>

int core = 0;
int direct = 0;


typedef struct thd
{
  int id;
  int core;
} thd_t;

void* 
test(void* d)
{
  thd_t* td = (thd_t*) d;
  int id = td->id;
  int my_core = the_cores[core + td->core];
  if (direct)
    {
      my_core = core + td->core;
    }

  printf("[%-2d] running on core %d\n", id, my_core);
  set_cpu(my_core);

  volatile long long unsigned int i, j, k = 0;
  while (1)
    {
      for (i = 2; i < 2e32; i++)
	{
	  for (j = 2; j < i; j++)
	    {
	      if (i % j == 0)
		{
		  k--;
		  break;
		}
	    }
	  k++;
	}
    }
}


int
main(int a, char** v)
{
  printf("//Usage: ./noise [NUM_THREADS] [THREADS_PER_CORE] [DIRECT]\n"
	 "//        where DIRECT selects whether to start with set_cpu(0), or set_cpu(the_cores[0])\n");
  int num_threads = 1;
  int per_core = 1;

  if (a > 1)
    {
      num_threads = atoi(v[1]);
    }
  if (a > 2)
    {
      per_core = atoi(v[2]);
    }
  if (a > 3)
    {
      direct = 1;
    }

  printf("# num of threads: %d / threads per core: %d / direct placement: %d\n", num_threads, per_core, direct);


  pthread_t threads[num_threads];
  pthread_attr_t attr;
  int rc;
  void *status;
    
  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
  thd_t* d = (thd_t*) malloc(num_threads * sizeof(thd_t));
  assert(d != NULL);

  int c = 0;
  int this = 0;

  long t;
  for(t = 0; t < num_threads; t++)
    {
      d[t].id = t;
      d[t].core = c;
      rc = pthread_create(&threads[t], &attr, test, d + t);
      if (rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}

      if (++this == per_core)
	{
	  this = 0;
	  c++;
	}
   
    }
    

  /* Free attribute and wait for the other threads */
  pthread_attr_destroy(&attr);
    
  for(t = 0; t < num_threads; t++) 
    {
      rc = pthread_join(threads[t], &status);
      if (rc) 
	{
	  printf("ERROR; return code from pthread_join() is %d\n", rc);
	  exit(-1);
	}
    }

  free(d);

  return 0;
}
