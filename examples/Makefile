all: arc_cavity ellipse_cavity trapping wedge triangle disknormcond kitenormcond invellipse disk

trapping: trapping.o
	$(CPP) $(CFLAGS) -o ../bin/examples/trapping trapping.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

wedge: wedge.o
	$(CPP) $(CFLAGS) -o ../bin/examples/wedge wedge.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

invellipse: invellipse.o
	$(CPP) $(CFLAGS) -o ../bin/examples/invellipse invellipse.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)


triangle: triangle.o
	$(CPP) $(CFLAGS) -o ../bin/examples/triangle triangle.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

disknormcond: disknormcond.o
	$(CPP) $(CFLAGS) -o ../bin/examples/disknormcond disknormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

kitenormcond: kitenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/examples/kitenormcond kitenormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

disk: disk.o
	$(CPP) $(CFLAGS) -o ../bin/examples/disk disk.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

ellipse_cavity: ellipse_cavity.o
	$(CPP) $(CFLAGS) -o ../bin/examples/ellipse_cavity ./ellipse_cavity.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

arc_cavity: arc_cavity.o
	$(CPP) $(CFLAGS) -o ../bin/examples/arc_cavity ./arc_cavity.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

