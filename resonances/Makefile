all: trapresonance

trapresonance: trapresonance.o
	$(CPP) $(CFLAGS) -o ../bin/resonances/trapresonance trapresonance.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)


%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

