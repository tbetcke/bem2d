include ../make.inc

OBJECTS	= Point.o

all: $(OBJECTS)
	ar rcs bem2d.a $(OBJECTS)

###	$(CPP) $(CFLAGS) -o $(NAME) $(OBJECTS) $(LDFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) -I$(GSLINC) -I$(ATLASINC) -c $<

