SMC_ROOT = /nethome/rvuduc3/local/jinx/bebop-smc

CC = icc
CFLAGS = -std=c99 -I$(SMC_ROOT)/include -I.
COPTFLAGS = -O3 -g
LDFLAGS = -Wl,-rpath -Wl,$(SMC_ROOT)/lib -L$(SMC_ROOT)/lib -lsparse_matrix_converter -lbebop_util
CLEANFILES =

#=======================================================================

all: cg$(EXEEXT)

#=======================================================================

cg_HDRS = csr.h cg.h timer.h cycle.h scan.h
cg_SRCS = driver.c csr.c csr_matvec.c cg.c timer.c scan.c
cg_OBJS = $(cg_SRCS:.c=.o)

cg$(EXEEXT): $(cg_OBJS)
	$(CC) $(COPTFLAGS) -o $@ $(cg_OBJS) $(LDFLAGS)

CLEANFILES += cg$(EXEEXT) $(cg_OBJS)
CLEANFILES += rhist__sequential.out
CLEANFILES += rhist__parfor.out
CLEANFILES += rhist__recspawn.out
CLEANFILES += rhist__segscan.out

#=======================================================================

%.o: %.c
	$(CC) $(CFLAGS) $(COPTFLAGS) -o $@ -c $<

#=======================================================================

clean:
	rm -rf core *~ $(CLEANFILES)

# eof
