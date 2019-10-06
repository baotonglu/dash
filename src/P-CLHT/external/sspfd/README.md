sspfd
=====

sspfd if a super-simple profiler, able to measure timings with granularity of a single CPU cycle. sspfd gathers the timings and when the application decides to print the data, it does some statistical analysis on the results (i.e., std deviation, abs deviation, clustering of values) and prints the results.  
For example, you can use sspfd to measure the latency of acquiring a lock, or of performing an atomic operation on some data.

* Website             : https://github.com/trigonak/sspfd.git
* Author              : Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>

Installation:
-------------

Compile the library using `make` in the base folder.

Test the installation using `sspfd_test` (see the available options with `./sspfd_test -h`).

Use the interface in `sspfd.h` and link your application with the `libsspfd.a` library.

Interface:
----------

`SSPFD_DO_TIMINGS` should be set to 1 if you want to take measurements, else you can set it to any other value to remove `sspfd`'s functionality (and overhead).

For an example of how to use `sspfd`, refer to `sspfd_test.c`.

Interface of sspfd:

* `SSPFDINIT(num_stores, num_entries, id)` : initialize `num_stores` stores, with `num_entries` number of measurement entries each. Use `id` as the thread id for printing purposes. Each store should be used to measure a single event.
* `SSPFDTERM()` : terminate (free) the initialized stores.
* `SSPFDI(store)` : start measurement (i.e., take start timestamp at this point) for store `store`.
* `SSPFDO(store, entry)` : stop measuring (i.e., take stop timestamp at this point) for store `store` and store the duration since `SSPFDI(store)` in entry `entry`. 
* `SSPFDSTATS(store, num_ops, statsp)` : generate statistics for the measurements in store `store` for the first `num_ops` values. Store the results in `statsp` pointer to a `sspfd_stats_t` structure. 
* `SSPFDPRINT(statsp)` : print the statistics in `statsp` pointer to a `sspfd_stats_t` structure. 
* `SSPFDPRINTV(store, num_print)` : print the first `num_print` measurements from store `store`.
* `SSPFDP(store, num_vals)` : generate statistics and print them for the first `num_vals` values of store `store`.
* `SSPFDPN(store, num_vals, num_print)` : generate statistics and print them for the first `num_vals` values of store `store`. Additionally, print the first `num_print` measurements of this store.
* `SSPFDPREFTCH(store, entry)` : prefetch entry `entry` for store `store`, so that the overheads are minimized (only necessary if the application has accessed a huge amount of data before using `sspfd`.

Interpreting the results:
-------------------------

The comments prefixed with "#######" explain the results.

<pre>
####### warning for indicating that the correction for sspfd was not stable
* warning: avg sspfd correction is 64.6 with std deviation: 47.4%. Recalculating.

####### correction = the cost of a getticks (to get a timestamp) + the cost of storing
####### the timestamps and calculating the interval between two timestamps
 -- sspfd correction: 64 (std deviation: 1.6%)

####### results
 ---- statistics:
####### global avg, abs & std dev, num of measurements
[00]     avg : 1.7        abs dev : 0.8        std dev : 30.4       num     : 100000
####### min & max values
[00]     min : 0.0        (element:      2)    max     : 9614.0     (element:  92787)
####### clustering of values around the global avg
####### 0-10% contains values that are in the range [(0.9*avg)..(1.1*avg)]
[00]   0-10% : 0          (  0.0%  |  avg:    -nan  |  abs dev:   -nan  |  std dev:   -nan =  -nan%)
[00]  10-25% : 38214      ( 38.2%  |  avg:     2.0  |  abs dev:    0.0  |  std dev:    0.0 =   0.0%)
[00]  25-50% : 43422      ( 43.4%  |  avg:     1.0  |  abs dev:    0.0  |  std dev:    0.0 =   0.0%)
[00]  50-75% : 0          (  0.0%  |  avg:    -nan  |  abs dev:   -nan  |  std dev:   -nan =  -nan%)
[00] 75-100% : 18364      ( 18.4%  |  avg:     2.6  |  abs dev:    1.7  |  std dev:   71.0 = 2741.5%)
</pre>
