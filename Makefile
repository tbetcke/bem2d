include make.inc

ifeq ($(debug),y)
   MAKE += debug=y
endif

all:	bem2d tests


.PHONY: tests
tests: bem2d
	cd tests; $(MAKE); ./testall

.PHONY: bem2d
bem2d:
	cd libs; $(MAKE)
 
clean:
	rm -f libs/*.o; rm -f tests/*.o

