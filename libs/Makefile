include ../config/$(config).inc

ifeq ($(debug),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif


OBJECTS	= Point.o Element.o Basis.o quadrature.o QuadPoints.o

all: $(OBJECTS)
	ar rcs bem2d.a $(OBJECTS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -I$(GSLINC) -I$(ATLASINC) -c $<

clean:
	rm -f *.o *~ bem2d.a
