ifdef config
  include ./config/$(config).inc
else
  include ./config/default.inc
endif

ifeq ($(debug),y)
 CFLAGS += $(DEBUGFLAGS)
else
 CFLAGS += $(OPTFLAGS)
endif

ifeq ($(profile),y)
   CFLAGS += -pg
endif

ifeq ($(mpi),y)
   CFLAGS += -DBEM2DMPI
endif

ifeq ($(omp),y)
   CFLAGS += -fopenmp
endif

export CFLAGS
export LDFLAGS
export CPP

all:	bem2d tests examples


.PHONY: tests
tests: bem2d dirs
	cd tests; make

.PHONY: bem2d
bem2d:
	cd lib; make

.PHONY: clean 
clean:
	cd lib; make clean
	cd tests; make clean
	cd examples; make clean
	cd normcomps; make clean
	rm -f *~ *.o
	rm -rf bin

.PHONY: examples
examples: bem2d dirs
	cd examples; make

.PHONY: dirs
dirs:
	test -d bin || mkdir bin;
	test -d bin/examples || mkdir bin/examples;

.PHONY: normcomps
normcomps: bem2d dirs
	test -d bin/normcomps || mkdir bin/normcomps
	cd normcomps; make

.PHONY: numrangecomps
numrangecomps: bem2d dirs
	test -d bin/numrangecomps || mkdir bin/numrangecomps
	cd numrangecomps; make


.PHONY: commit
commit: clean
	git commit --all
