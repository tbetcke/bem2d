include ../config/$(config).inc

ifeq ($(debug),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif

ifeq ($(profile),y)
   MAKE += profile=y
   CFLAGS += -pg
endif


all: $(OBJECTS)
	$(CPP) $(CFLAGS) -o ../bin/test test.cpp ../lib/bem2d.a $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~
