include ../make.inc

ifeq ($(debug),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif


OBJECTS	= Point.o Element.o

all: $(OBJECTS)
	ar rcs bem2d.a $(OBJECTS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -I$(GSLINC) -I$(ATLASINC) -c $<

clean:
	rm -f *.o *~ bem2d.a
