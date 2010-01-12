all: disknormcond kitenormcond ellipsenormcond squarenormcond cracknormcond rectanglenormcond trappingnormcond

disknormcond: disknormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/disknormcond disknormcond.cpp ../lib/bem2d.a $(LDFLAGS)

kitenormcond: kitenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/kitenormcond kitenormcond.cpp ../lib/bem2d.a $(LDFLAGS)

ellipsenormcond: ellipsenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/ellipsenormcond ellipsenormcond.cpp ../lib/bem2d.a $(LDFLAGS)

squarenormcond: squarenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/squarenormcond squarenormcond.cpp ../lib/bem2d.a $(LDFLAGS)

cracknormcond: cracknormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/cracknormcond cracknormcond.cpp ../lib/bem2d.a $(LDFLAGS)

rectanglenormcond: rectanglenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/rectanglenormcond rectanglenormcond.cpp ../lib/bem2d.a $(LDFLAGS)

trappingnormcond: trappingnormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/trappingnormcond trappingnormcond.cpp ../lib/bem2d.a $(LDFLAGS)



%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

