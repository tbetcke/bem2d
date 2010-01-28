all: disk trapping

disk: disk.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/disk disk.cpp ../lib/bem2d.a $(LDFLAGS)

trapping: trapping.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/trapping trapping.cpp ../lib/bem2d.a $(LDFLAGS)


%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

