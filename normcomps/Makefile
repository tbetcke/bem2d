all: disknormcond kitenormcond ellipsenormcond squarenormcond cracknormcond rectanglenormcond trappingnormcond ellipse_cavity_normcond arc_cavity_normcond

disknormcond: disknormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/disknormcond disknormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

kitenormcond: kitenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/kitenormcond kitenormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

ellipsenormcond: ellipsenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/ellipsenormcond ellipsenormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

squarenormcond: squarenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/squarenormcond squarenormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

cracknormcond: cracknormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/cracknormcond cracknormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

rectanglenormcond: rectanglenormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/rectanglenormcond rectanglenormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

trappingnormcond: trappingnormcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/trappingnormcond trappingnormcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

ellipse_cavity_normcond: ellipse_cavity_normcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/ellipse_cavity_normcond ellipse_cavity_normcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)

arc_cavity_normcond: arc_cavity_normcond.o
	$(CPP) $(CFLAGS) -o ../bin/normcomps/arc_cavity_normcond arc_cavity_normcond.cpp ../lib/bem2d.a ../amos/amos.a $(LDFLAGS)


%.o: %.cpp
	$(CPP) $(CFLAGS)  -c $<

clean: 
	rm -f *.o *~ *.s *.ii

