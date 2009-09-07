include ./config/$(config).inc
MAKE += config=$(config)

ifeq ($(debug),y)
   MAKE += debug=y
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif


all:	bem2d tests


.PHONY: tests
tests: bem2d
	cd tests; $(MAKE)

.PHONY: bem2d
bem2d:
	cd libs; $(MAKE)

.PHONY: clean 
clean:
	cd libs; $(MAKE) clean
	cd tests; $(MAKE) clean
	rm -f *~ *.o test

testfile: bem2d tests
	$(CPP) $(CFLAGS) -I./libs -o test test.cpp ./libs/bem2d.a $(LDFLAGS)


