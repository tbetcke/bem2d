all: $(OBJECTS)
	$(CPP) $(CFLAGS) -o ../bin/test test.cpp ../lib/bem2d.a ../amos/amos.a x$(LDFLAGS) 
%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii
