OBJECTS	= bem2d_point.o bem2d_element.o bem2d_basis.o bem2d_quadrature.o bem2d_quadpoints.o bem2d_kernel.o bem2d_geometry.o \
		  bem2d_soundsoftscattering.o bem2d_mathroutines.o bem2d_outputroutines.o bem2d_outputhandler.o bem2d_fun.o \
		  bem2d_shape.o bem2d_curve.o bem2d_basistypes.o bem2d_pnpoly.o bem2d_mpi.o bem2d_curvetypes.o

all: $(OBJECTS)
	ar rcs bem2d.a $(OBJECTS)

%.o: %.f
	$(F77) $(FFLAGS) -c $<

%.o: %.cpp
	$(CPP) $(CFLAGS) -c $<

clean:
	rm -f *.o *~ *.s *.ii bem2d.a
