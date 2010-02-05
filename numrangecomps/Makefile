all: disk diskhconv trapping invellipse lshape lshapehconv squarehconv

disk: disk.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/disk disk.cpp ../lib/bem2d.a $(LDFLAGS)

diskhconv: diskhconv.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/diskhconv diskhconv.cpp ../lib/bem2d.a $(LDFLAGS)


trapping: trapping.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/trapping trapping.cpp ../lib/bem2d.a $(LDFLAGS)

invellipse: invellipse.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/invellipse invellipse.cpp ../lib/bem2d.a $(LDFLAGS)

lshape: lshape.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/lshape lshape.cpp ../lib/bem2d.a $(LDFLAGS)

lshapehconv: lshapehconv.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/lshapehconv lshapehconv.cpp ../lib/bem2d.a $(LDFLAGS)

squarehconv: squarehconv.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/squarehconv squarehconv.cpp ../lib/bem2d.a $(LDFLAGS)


%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

