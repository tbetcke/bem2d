all: $(OBJECTS)
	$(CPP) $(CFLAGS) -o ../bin/test test.cpp ../lib/bem2d.a $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii
