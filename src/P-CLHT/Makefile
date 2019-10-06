#################################
# Architecture dependent settings
#################################

CFLAGS = -D_GNU_SOURCE -DCLFLUSH

ifeq ($(DEBUG),1)
  DEBUG_FLAGS=-Wall -ggdb -g -DDEBUG
  CFLAGS += -O0 -DADD_PADDING -fno-inline
else ifeq ($(DEBUG),2)
  DEBUG_FLAGS=-Wall
  CFLAGS += -O0 -DADD_PADDING -fno-inline
else ifeq ($(DEBUG),3)
  DEBUG_FLAGS=-Wall -g -ggdb 
  CFLAGS += -O3 -DADD_PADDING -fno-inline
else
  DEBUG_FLAGS=-Wall
  CFLAGS += -O3 -DADD_PADDING
endif

ifeq ($(SET_CPU),0)
	CFLAGS += -DNO_SET_CPU
endif

ifeq ($(LATENCY),1)
	CFLAGS += -DCOMPUTE_LATENCY -DDO_TIMINGS
endif

ifeq ($(LATENCY),2)
	CFLAGS += -DCOMPUTE_LATENCY -DDO_TIMINGS -DUSE_SSPFD -DLATENCY_ALL_CORES=0
	LIBS += $(SSPFD) -lm
endif

ifeq ($(LATENCY),3)
	CFLAGS += -DCOMPUTE_LATENCY -DDO_TIMINGS -DUSE_SSPFD -DLATENCY_ALL_CORES=1
	LIBS += $(SSPFD) -lm
endif

ifeq ($(LATENCY),4)
	CFLAGS += -DCOMPUTE_LATENCY -DDO_TIMINGS -DUSE_SSPFD -DLATENCY_PARSING=1
	LIBS += $(SSPFD) -lm
endif

ifeq ($(LATENCY),5)
	CFLAGS += -DCOMPUTE_LATENCY -DDO_TIMINGS -DUSE_SSPFD -DLATENCY_PARSING=1 -DLATENCY_ALL_CORES=1
	LIBS += $(SSPFD) -lm
endif

TOP := $(patsubst %/,%,$(dir $(lastword $(MAKEFILE_LIST))))

LIBS+=-L$(TOP)/external/lib -L$(TOP)

SRCPATH := $(TOP)/src
MAININCLUDE := $(TOP)/include

ifeq ($(M),1)
LIBS += -lsspfd
CFLAGS += -DUSE_SSPFD
endif

# default setings
PLATFORM=-DDEFAULT #-DPERFCNT
GCC=gcc
PLATFORM_NUMA=0
OPTIMIZE=
LIBS += -lrt -lpthread -lm  -lclht -lssmem

UNAME := $(shell uname -n)

ifeq ($(UNAME), lpd48core)
PLATFORM=-DOPTERON
GCC=gcc-4.8
PLATFORM_NUMA=1
OPTIMIZE=-DOPTERON_OPTIMIZE
PLATFORM_NUMA=1
endif

ifeq ($(UNAME), lpdxeon2680)
PLATFORM=-DXEON2
GCC=gcc
PLATFORM_NUMA=1
OPTIMIZE=
endif

ifeq ($(UNAME), lpdpc4)
PLATFORM=-DCOREi7
GCC=gcc
PLATFORM_NUMA=0
OPTIMIZE=
endif

ifeq ($(UNAME), lpdpc34)
PLATFORM=-DCOREi7 -DRTM
GCC=gcc-4.8
PLATFORM_NUMA=0
OPTIMIZE=
LIBS += -mrtm
endif

ifeq ($(UNAME), diascld9)
PLATFORM=-DOPTERON2
GCC=gcc
LIBS += 
endif

ifeq ($(UNAME), diassrv8)
PLATFORM=-DXEON
GCC=gcc
PLATFORM_NUMA=1
LIBS += -lnuma
endif

ifeq ($(UNAME), diascld19)
PLATFORM=-DXEON2
GCC=gcc
LIBS += 
endif

ifeq ($(UNAME), maglite)
PLATFORM=-DSPARC
GCC:=/opt/csw/bin/gcc
LIBS+= 
CFLAGS += -m64 -mcpu=v9 -mtune=v9
endif

ifeq ($(UNAME), parsasrv1.epfl.ch)
PLATFORM=-DTILERA
GCC=tile-gcc
LIBS += -ltmc
endif

ifeq ($(UNAME), smal1.sics.se)
PLATFORM=-DTILERA
GCC=tile-gcc
LIBS += -ltmc
endif

ifeq ($(UNAME), ol-collab1)
PLATFORM=-DT44
GCC=/usr/sfw/bin/gcc
CFLAGS += -m64
LIBS+= 
endif

ifeq ($(PLATFORM_NUMA), 1)
CFLAGS += -DNUMA
LIBS += -lnuma
endif

CFLAGS += $(PLATFORM)
CFLAGS += $(OPTIMIZE)
CFLAGS += $(DEBUG_FLAGS)

INCLUDES := -I$(MAININCLUDE) -I$(TOP)/external/include
OBJ_FILES := clht_gc.o

SRC := src

BMARKS := bmarks

MAIN_BMARK := $(BMARKS)/test.c     # no memory allocation
#MAIN_BMARK := $(BMARKS)/test_ro.c  # read-only benchmark
#MAIN_BMARK := $(BMARKS)/test_mem.c  # memory allocation

ALL = 	clht_lb clht_lb_res clht_lb_res_no_next clht_lb_ro clht_lb_linked clht_lb_packed \
	clht_lf clht_lf_res clht_lf_only_map_rem

REST = 	math_cache_lb math_cache_lf \
	snap_stress noise

default: normal

dependencies:
	./scripts/make_dependencies.sh

rest: $(REST)

all: $(ALL)

.PHONY: $(ALL) \
	libclht_lb.a libclht_lb_res.a libclht_lb_res_no_next.a \
	libclht_lb_linked.a libclht_lb_packed.a libclht_lb_lock_ins.a \
	libclht_lf.a libclht_lf_res.a libclht_lf_only_map_rem.a

normal: clht_lb_res clht_lf_res 


%.o:: $(SRC)/%.c 
	$(GCC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

clht_gc_linked.o: $(SRC)/clht_gc.c
	$(GCC) -DCLHT_LINKED $(CFLAGS) $(INCLUDES) -o clht_gc_linked.o -c $(SRC)/clht_gc.c

################################################################################
# library
################################################################################

TYPE = clht_lb
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lb.o $(OBJ_FILES)

TYPE = clht_lb_res
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lb_res.o $(OBJ_FILES)

TYPE = clht_lb_res_no_next
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lb_res_no_next.o $(OBJ_FILES)

TYPE = clht_lb_linked
OBJ = $(TYPE).o
lib$(TYPE).a: clht_gc_linked.o $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lb_linked.o clht_gc_linked.o

TYPE = clht_lb_packed
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lb_packed.o $(OBJ_FILES)

TYPE = clht_lb_lock_ins
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lb_lock_ins.o $(OBJ_FILES)

TYPE = clht_lf
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lf.o $(OBJ_FILES)

TYPE = clht_lf_res
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lf_res.o $(OBJ_FILES)

TYPE = clht_lf_only_map_rem
OBJ = $(TYPE).o
lib$(TYPE).a: $(OBJ_FILES) $(OBJ)
	@echo Archive name = libclht.a
	ar -d libclht.a *
	ar -r libclht.a clht_lf_only_map_rem.o $(OBJ_FILES)

################################################################################
# lock-based targets
################################################################################

TYPE = clht_lb
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) -DNO_RESIZE $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) -o clht_lb $(LIBS) -lclht

TYPE = clht_lb_res
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) $(INCLUDES) $(CFLAGS) $(MAIN_BMARK) $(SRC)/clht_lb_res.c -o clht_lb_res $(LIBS)

TYPE = clht_lb_res_no_next
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) -o clht_lb_nn $(LIBS)

TYPE = clht_lb_linked
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) -o clht_lb_linked $(LIBS)

TYPE = clht_lb_packed
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) $(SRC)/clht_lb_packed.c -o clht_lb_packed $(LIBS)

TYPE = clht_lb_lock_ins
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) -DLOCK_INS $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) -o clht_lb_lock_ins $(LIBS)

################################################################################
# lock-free targets
################################################################################

TYPE = clht_lf
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) -DLOCKFREE $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) -o clht_lf $(LIBS)

TYPE = clht_lf_res
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) -DLOCKFREE_RES $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) -o clht_lf_res $(LIBS)

TYPE = clht_lf_only_map_rem
$(TYPE): $(MAIN_BMARK) lib$(TYPE).a
	$(GCC) -DLOCKFREE $(CFLAGS) $(INCLUDES) $(MAIN_BMARK) -o clht_lf_only_map_rem $(LIBS)

################################################################################
# other tests
################################################################################

math_cache_lb: $(BMARKS)/math_cache.c libclht_lb_res.a
	$(GCC) $(CFLAGS) $(INCLUDES) $(BMARKS)/math_cache.c -o math_cache_lb $(LIBS)

math_cache_lf: $(BMARKS)/math_cache.c libclht_lf_res.a
	$(GCC) -DLOCKFREE $(CFLAGS) $(INCLUDES) $(BMARKS)/math_cache.c -o math_cache_lf $(LIBS)

snap_stress: $(BMARKS)/snap_stress.c libclht_lf_res.a
	$(GCC) -DLOCKFREE $(CFLAGS) $(INCLUDES) $(BMARKS)/snap_stress.c -o snap_stress $(LIBS)

noise: $(BMARKS)/noise.c $(OBJ_FILES)
	$(GCC) $(CFLAGS) $(INCLUDES) $(OBJ_FILES) $(BMARKS)/noise.c -o noise $(LIBS)



clean:				
	rm -f *.o *.a clht_* math_cache* snap_stress
