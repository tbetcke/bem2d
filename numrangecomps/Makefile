all: disk diskhconv trapping invellipse lshape lshapehconv squarehconv square triangle squaredensity gww poly kite ellipse trappingresrange

disk: disk.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/disk disk.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

diskhconv: diskhconv.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/diskhconv diskhconv.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)


trapping: trapping.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/trapping trapping.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

invellipse: invellipse.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/invellipse invellipse.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

lshape: lshape.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/lshape lshape.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

lshapehconv: lshapehconv.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/lshapehconv lshapehconv.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

squarehconv: squarehconv.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/squarehconv squarehconv.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

square: square.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/square square.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

triangle: triangle.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/triangle triangle.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

squaredensity: squaredensity.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/squaredensity squaredensity.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

gww: gww.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/gww gww.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

kite: kite.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/kite kite.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

poly: poly.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/poly poly.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

ellipse: ellipse.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/ellipse ellipse.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

trappingresrange: trappingresrange.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/trappingresrange trappingresrange.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)


%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

