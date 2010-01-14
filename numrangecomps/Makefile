all: disknormcond

disknormcond: disk.o
	$(CPP) $(CFLAGS) -o ../bin/numrangecomps/disk disk.cpp ../lib/bem2d.a $(LDFLAGS)

trappingnormcond: trapping.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/trapping trapping.cpp ../lib/bem2d.a $(LDFLAGS)


%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

