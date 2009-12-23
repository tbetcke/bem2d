all: trapping wedge triangle

trapping: trapping.o
	$(CPP) $(CFLAGS) -o ../bin/examples/trapping trapping.cpp ../lib/bem2d.a $(LDFLAGS)

wedge: wedge.o
	$(CPP) $(CFLAGS) -o ../bin/examples/wedge wedge.cpp ../lib/bem2d.a $(LDFLAGS)

triangle: triangle.o
	$(CPP) $(CFLAGS) -o ../bin/examples/triangle triangle.cpp ../lib/bem2d.a $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

