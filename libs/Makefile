include ../config/$(config).inc

ifeq ($(debug),y)
   CFLAGS += $(DEBUGFLAGS)
else
   CFLAGS += $(OPTFLAGS)
endif


OBJECTS	= bem2d_point.o bem2d_element.o bem2d_basis.o bem2d_quadrature.o bem2d_quadpoints.o bem2d_kernel.o bem2d_geometry.o \
		  bem2d_soundsoftscattering.o bem2d_mathroutines.o bem2d_outputroutines.o bem2d_outputhandler.o bem2d_fun.o \
		  bem2d_shape.o bem2d_curve.o bem2d_basistypes.o

all: $(OBJECTS)
	ar rcs bem2d.a $(OBJECTS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -I$(GSLINC) -I$(ATLASINC) -c $<

clean:
	rm -f *.o *~ bem2d.a
