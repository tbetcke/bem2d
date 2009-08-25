include ../make.inc

ifeq ($(debug),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif

LDFLAGS += -lboost_unit_test_framework

OBJECTS	= runner.o testPoint.o

all: $(OBJECTS)
	$(CPP) $(CFLAGS) -o testall $(OBJECTS) ../libs/bem2d.a $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -I$(GSLINC) -I$(ATLASINC) -c $<

clean: 
	rm -f *.o
