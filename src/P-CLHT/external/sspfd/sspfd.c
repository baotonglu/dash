/*   
 *   File: sspfd.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: a fine-grained profiler based on rdtsc
 *   sspfd.c is part of sspfd
 *
 * The MIT License (MIT)
 *
 * Copyright (C) 2013  Vasileios Trigonakis
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

#include "sspfd.h"

__thread volatile size_t sspfd_num_stores;
__thread volatile ticks** sspfd_store;
__thread volatile ticks* _sspfd_s;
__thread ticks _sspfd_s_global;
__thread volatile ticks sspfd_correction;
__thread size_t SSPFD_ID = 0;

void 
sspfd_store_init(size_t num_stores, size_t num_entries, size_t id)
{
  SSPFD_ID = id;
  sspfd_num_stores = num_stores;

  _sspfd_s = (volatile ticks*) malloc(num_stores * sizeof(ticks));
  sspfd_store = (volatile ticks**) malloc(num_stores * sizeof(ticks*));
  assert(_sspfd_s != NULL && sspfd_store != NULL);

  volatile uint32_t i;
  for (i = 0; i < num_stores; i++)
    {
      sspfd_store[i] = (ticks*) malloc(num_entries * sizeof(ticks));
      assert(sspfd_store[i] != NULL);
      PREFETCHW((void*) &sspfd_store[i][0]);
    }

  int32_t tries = 10;
  uint32_t print_warning = 0;


  /* "enforcing" max freq if freq scaling is enabled */
  volatile uint64_t speed;
  for (speed = 0; speed < 20e6; speed += 4)
    {
      speed--;
    }

  sspfd_correction = 0;

#define SSPFD_CORRECTION_CONF 3
 retry:
  for (i = 0; i < num_entries; i++)
    {
      SSPFDI(0);
      asm volatile ("");
      SSPFDO(0, i);
    }

  sspfd_stats_t ad;
  sspfd_get_stats(0, num_entries, &ad);
  double std_pp = 100 * (1 - (ad.avg - ad.std_dev) / ad.avg);

  if (std_pp > SSPFD_CORRECTION_CONF)
    {
      if (print_warning++ == 1)	/* print warning if 2 failed attempts */
	{
#if defined(DEBUG)
	  printf("* warning: avg sspfd correction is %.1f with std deviation: %.1f%%. Recalculating.\n", 
		 ad.avg, std_pp);
#endif
	}
      if (tries-- > 0)
	{
	  goto retry;
	}
      else
	{
	  printf("* warning: setting sspfd correction with a high std deviation\n");
	}
    }

  sspfd_correction = ad.avg;
  assert(sspfd_correction > 0);
  
#if defined(DEBUG)
  printf(" -- sspfd correction: %llu (std deviation: %.1f%%)\n", (long long unsigned int) sspfd_correction, std_pp);
#endif
}

void
sspfd_store_term()
{
  volatile uint32_t i;
  for (i = 0; i < sspfd_num_stores; i++)
    {
      free((void*) sspfd_store[i]);
    }

  free((void*) _sspfd_s);
  free(sspfd_store);
}

inline void
sspfd_set_id(size_t id)
{
  SSPFD_ID = id;
}

inline size_t
sspfd_get_id()
{
  return SSPFD_ID;
}

static inline double
sspfd_absd(double x)
{
  if (x >= 0)
    {
      return x;
    }
  else 
    {
      return -x;
    }
}

#define llu long long unsigned int
void 
sspfd_print_stats(const sspfd_stats_t* sspfd_stats)
{
  printf("\n ---- statistics:\n");
  SSPFD_PRINT("    avg : %-10.1f abs dev : %-10.1f std dev : %-10.1f num     : %llu", 
	sspfd_stats->avg, sspfd_stats->abs_dev, sspfd_stats->std_dev, (llu) sspfd_stats->num_vals);
  SSPFD_PRINT("    min : %-10.1f (element: %6llu)    max     : %-10.1f (element: %6llu)", sspfd_stats->min_val, 
	(llu) sspfd_stats->min_val_idx, sspfd_stats->max_val, (llu) sspfd_stats->max_val_idx);
  double v10p = 100 * 
    (1 - (sspfd_stats->num_vals - sspfd_stats->num_dev_10p) / (double) sspfd_stats->num_vals);
  double std_10pp = 100 * (1 - (sspfd_stats->avg_10p - sspfd_stats->std_dev_10p) / sspfd_stats->avg_10p);
  SSPFD_PRINT("  0-10%% : %-10u (%5.1f%%  |  avg:  %6.1f  |  abs dev: %6.1f  |  std dev: %6.1f = %5.1f%%)", 
	sspfd_stats->num_dev_10p, v10p, sspfd_stats->avg_10p, sspfd_stats->abs_dev_10p, sspfd_stats->std_dev_10p, std_10pp);
  double v25p = 100 
    * (1 - (sspfd_stats->num_vals - sspfd_stats->num_dev_25p) / (double) sspfd_stats->num_vals);
  double std_25pp = 100 * (1 - (sspfd_stats->avg_25p - sspfd_stats->std_dev_25p) / sspfd_stats->avg_25p);
  SSPFD_PRINT(" 10-25%% : %-10u (%5.1f%%  |  avg:  %6.1f  |  abs dev: %6.1f  |  std dev: %6.1f = %5.1f%%)", 
	sspfd_stats->num_dev_25p, v25p, sspfd_stats->avg_25p, sspfd_stats->abs_dev_25p, sspfd_stats->std_dev_25p, std_25pp);
  double v50p = 100 * 
    (1 - (sspfd_stats->num_vals - sspfd_stats->num_dev_50p) / (double) sspfd_stats->num_vals);
  double std_50pp = 100 * (1 - (sspfd_stats->avg_50p - sspfd_stats->std_dev_50p) / sspfd_stats->avg_50p);
  SSPFD_PRINT(" 25-50%% : %-10u (%5.1f%%  |  avg:  %6.1f  |  abs dev: %6.1f  |  std dev: %6.1f = %5.1f%%)", 
	sspfd_stats->num_dev_50p, v50p, sspfd_stats->avg_50p, sspfd_stats->abs_dev_50p, sspfd_stats->std_dev_50p, std_50pp);
  double v75p = 100 * 
    (1 - (sspfd_stats->num_vals - sspfd_stats->num_dev_75p) / (double) sspfd_stats->num_vals);
  double std_75pp = 100 * (1 - (sspfd_stats->avg_75p - sspfd_stats->std_dev_75p) / sspfd_stats->avg_75p);
  SSPFD_PRINT(" 50-75%% : %-10u (%5.1f%%  |  avg:  %6.1f  |  abs dev: %6.1f  |  std dev: %6.1f = %5.1f%%)", 
	sspfd_stats->num_dev_75p, v75p, sspfd_stats->avg_75p, sspfd_stats->abs_dev_75p, sspfd_stats->std_dev_75p, std_75pp);
  double vrest = 100 * 
    (1 - (sspfd_stats->num_vals - sspfd_stats->num_dev_rst) / (double) sspfd_stats->num_vals);
  double std_rspp = 100 * (1 - (sspfd_stats->avg_rst - sspfd_stats->std_dev_rst) / sspfd_stats->avg_rst);
  SSPFD_PRINT("75-100%% : %-10u (%5.1f%%  |  avg:  %6.1f  |  abs dev: %6.1f  |  std dev: %6.1f = %5.1f%%)\n", 
	sspfd_stats->num_dev_rst, vrest, sspfd_stats->avg_rst, sspfd_stats->abs_dev_rst, sspfd_stats->std_dev_rst, std_rspp);
}

void
sspfd_get_stats(const size_t store, const size_t num_vals, sspfd_stats_t* sspfd_stats)
{
  volatile ticks* vals = sspfd_store[store];

  sspfd_stats->num_vals = num_vals;
  ticks sum_vals = 0;
  uint32_t i;
  for (i = 0; i < num_vals; i++)
    {
      if ((int64_t) vals[i] < 0)
	{
	  vals[i] = 0;
	}
      sum_vals += vals[i];
    }

  double avg = sum_vals / (double) num_vals;
  sspfd_stats->avg = avg;
  double max_val = 0;
  double min_val = DBL_MAX;
  uint64_t max_val_idx = 0, min_val_idx = 0;
  uint32_t num_dev_10p = 0; ticks sum_vals_10p = 0; double dev_10p = 0.1 * avg;
  uint32_t num_dev_25p = 0; ticks sum_vals_25p = 0; double dev_25p = 0.25 * avg;
  uint32_t num_dev_50p = 0; ticks sum_vals_50p = 0; double dev_50p = 0.5 * avg;
  uint32_t num_dev_75p = 0; ticks sum_vals_75p = 0; double dev_75p = 0.75 * avg;
  uint32_t num_dev_rst = 0; ticks sum_vals_rst = 0;

  double sum_adev = 0;		/* abs deviation */
  double sum_stdev = 0;		/* std deviation */
  for (i = 0; i < num_vals; i++)
    {
      double diff = vals[i] - avg;
      double ad = sspfd_absd(diff);
      if (vals[i] > max_val)
	{
	  max_val = vals[i];
	  max_val_idx = i;
	}
      else if (vals[i] < min_val)
	{
	  min_val = vals[i];
	  min_val_idx = i;
	}

      if (ad <= dev_10p)
	{
	  num_dev_10p++;
	  sum_vals_10p += vals[i];
	}
      else if (ad <= dev_25p)
	{
	  num_dev_25p++;
	  sum_vals_25p += vals[i];
	}
      else if (ad <= dev_50p)
	{
	  num_dev_50p++;
	  sum_vals_50p += vals[i];
	}
      else if (ad <= dev_75p)
	{
	  num_dev_75p++;
	  sum_vals_75p += vals[i];
	}
      else
	{
	  num_dev_rst++;
	  sum_vals_rst += vals[i];
	}

      sum_adev += ad;
      sum_stdev += ad*ad;
    }
  sspfd_stats->min_val = min_val;
  sspfd_stats->min_val_idx = min_val_idx;
  sspfd_stats->max_val = max_val;
  sspfd_stats->max_val_idx = max_val_idx;
  sspfd_stats->num_dev_10p = num_dev_10p;
  sspfd_stats->num_dev_25p = num_dev_25p;
  sspfd_stats->num_dev_50p = num_dev_50p;
  sspfd_stats->num_dev_75p = num_dev_75p;
  sspfd_stats->num_dev_rst = num_dev_rst;

  sspfd_stats->avg_10p = sum_vals_10p / (double) num_dev_10p;
  sspfd_stats->avg_25p = sum_vals_25p / (double) num_dev_25p;
  sspfd_stats->avg_50p = sum_vals_50p / (double) num_dev_50p;
  sspfd_stats->avg_75p = sum_vals_75p / (double) num_dev_75p;
  sspfd_stats->avg_rst = sum_vals_rst / (double) num_dev_rst;

  double sum_adev_10p = 0, sum_adev_25p = 0, sum_adev_50p = 0, sum_adev_75p = 0, sum_adev_rst = 0;
  double sum_stdev_10p = 0, sum_stdev_25p = 0, sum_stdev_50p = 0, sum_stdev_75p = 0, sum_stdev_rst = 0;

  /* pass again to calculate the deviations for the 10/25..p */
  for (i = 0; i < num_vals; i++)
    {
      double diff = vals[i] - avg;
      double ad = sspfd_absd(diff);
      if (ad <= dev_10p)
	{
	  double diff = vals[i] - sspfd_stats->avg_10p;
	  double ad = sspfd_absd(diff);
	  sum_adev_10p += ad;
	  sum_stdev_10p += (ad*ad);
	}
      else if (ad <= dev_25p)
	{
	  double diff = vals[i] - sspfd_stats->avg_25p;
	  double ad = sspfd_absd(diff);
	  sum_adev_25p += ad;
	  sum_stdev_25p += (ad*ad);
	}
      else if (ad <= dev_50p)
	{
	  double diff = vals[i] - sspfd_stats->avg_50p;
	  double ad = sspfd_absd(diff);
	  sum_adev_50p += ad;
	  sum_stdev_50p += (ad*ad);
	}
      else if (ad <= dev_75p)
	{
	  double diff = vals[i] - sspfd_stats->avg_75p;
	  double ad = sspfd_absd(diff);
	  sum_adev_75p += ad;
	  sum_stdev_75p += (ad*ad);
	}
      else
	{
	  double diff = vals[i] - sspfd_stats->avg_rst;
	  double ad = sspfd_absd(diff);
	  sum_adev_rst += ad;
	  sum_stdev_rst += (ad*ad);
	}
    }

  sspfd_stats->abs_dev_10p = sum_adev_10p / num_dev_10p; 
  sspfd_stats->abs_dev_25p = sum_adev_25p / num_dev_25p; 
  sspfd_stats->abs_dev_50p = sum_adev_50p / num_dev_50p; 
  sspfd_stats->abs_dev_75p = sum_adev_75p / num_dev_75p; 
  sspfd_stats->abs_dev_rst = sum_adev_rst / num_dev_rst; 

  sspfd_stats->std_dev_10p = sqrt(sum_stdev_10p / num_dev_10p); 
  sspfd_stats->std_dev_25p = sqrt(sum_stdev_25p / num_dev_25p); 
  sspfd_stats->std_dev_50p = sqrt(sum_stdev_50p / num_dev_50p); 
  sspfd_stats->std_dev_75p = sqrt(sum_stdev_75p / num_dev_75p); 
  sspfd_stats->std_dev_rst = sqrt(sum_stdev_rst / num_dev_rst); 

  double adev = sum_adev / num_vals;
  sspfd_stats->abs_dev = adev;
  double abs_dev_per = 100 * (1 - (sspfd_stats->avg - sspfd_stats->abs_dev) / sspfd_stats->avg);
  sspfd_stats->abs_dev_perc = abs_dev_per;
  double stdev = sqrt(sum_stdev / num_vals);
  sspfd_stats->std_dev = stdev;
  double std_dev_per = 100 * (1 - (sspfd_stats->avg - sspfd_stats->std_dev) / sspfd_stats->avg);
  sspfd_stats->std_dev_perc = std_dev_per;
}
