include ../config/$(config).inc

ifeq ($(debug),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif

LDFLAGS += -lboost_unit_test_framework-mt

OBJECTS	= UnitTest.o testPoint.o testConstElement.o

all: $(OBJECTS)
	$(CPP) $(CFLAGS) -o UnitTest $(OBJECTS) ../libs/bem2d.a $(LDFLAGS)
	$(CPP) $(CFLAGS) -o test test.cpp ../libs/bem2d.a $(LDFLAGS)
	./UnitTest

%.o: %.cpp
	$(CPP) $(CFLAGS) -I$(GSLINC) -I$(ATLASINC) -c $<

clean: 
	rm -f *.o *~ UnitTest
