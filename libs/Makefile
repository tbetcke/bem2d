include ../config/$(config).inc

ifeq ($(debug),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif


OBJECTS	= Point.o Element.o Basis.o quadrature.o QuadPoints.o kernel.o geometry.o \
		  soundsoftscattering.o mathroutines.o outputroutines.o Outputhandler.o Bem2dfun.o \
		  Shape.o Curve.o

all: $(OBJECTS)
	ar rcs bem2d.a $(OBJECTS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -I$(GSLINC) -I$(ATLASINC) -c $<

clean:
	rm -f *.o *~ bem2d.a
