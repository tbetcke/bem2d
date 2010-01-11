all: disknormcond

disknormcond: disknormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/disknormcond disknormcond.cpp ../lib/bem2d.a $(LDFLAGS)

kitenormcond: kitenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/examples/kitenormcond kitenormcond.cpp ../lib/bem2d.a $(LDFLAGS)


%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

