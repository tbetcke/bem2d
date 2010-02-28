OBJECTS	= dgamln.o xerror.o zabs.o zacai.o zacon.o \
zairy.o zasyi.o zbesh.o zbinu.o zbknu.o zbuni.o zbunk.o zdiv.o zexp.o \
zkscl.o zlog.o zmlri.o zmlt.o zrati.o zs1s2.o zseri.o zshch.o zsqrt.o \
zuchk.o zunhj.o zuni1.o zuni2.o zunik.o zunk1.o zunk2.o zuoik.o zwrsk.o \
i1mach.o d1mach.o


all: $(OBJECTS)
	ar rcs amos.a $(OBJECTS)

%.o: %.f
	$(F77) $(FFLAGS) -c $<

clean:
	rm -f *.o *~ *.s *.ii amos.a
