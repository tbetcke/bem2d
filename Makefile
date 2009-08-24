include make.inc

ifeq ($(DEBUG),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif

all:	bem2d tests


.PHONY: tests
tests: bem2d
	cd $(TESTDIR); $(MAKE)

.PHONY: bem2d
bem2d:
	cd $(SRCDIR); $(MAKE)
 
clean:
	rm -f $(SRCDIR)/*.o; rm -f $(TESTDIR)/*.o

