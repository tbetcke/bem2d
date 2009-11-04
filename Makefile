include ./config/$(config).inc
MAKE += config=$(config)

ifeq ($(debug),y)
   MAKE += debug=y
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif

ifeq ($(profile),y)
   MAKE += profile=y
   CFLAGS += -pg
endif


all:	bem2d tests


.PHONY: tests
tests: bem2d
	cd tests; $(MAKE)

.PHONY: bem2d
bem2d:
	cd $(SRCDIR); $(MAKE)

.PHONY: clean 
clean:
	cd lib; $(MAKE) clean
	cd tests; $(MAKE) clean
	rm -f *~ *.o
	rm -rf ./bin/*

testfile: bem2d tests
	$(CPP) $(CFLAGS) -I./libs -o test test.cpp ./libs/bem2d.a $(LDFLAGS)


