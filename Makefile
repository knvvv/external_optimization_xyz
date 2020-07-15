CC = g++
LDIR = ./lib
IDIR = ./include
CFLAGS = -fPIC -I$(IDIR) -L$(LDIR)
LIBS = -lgsl -lgslcblas
OBJ = cfiles/mainrand.o cfiles/cycledef.o cfiles/main.o

mylib.so: $(OBJ)
	$(CC) -o $@ $^ -shared $(CFLAGS) $(LIBS)

cfiles/main.o: cfiles/main.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

cfiles/cycledef.o: cfiles/cycledef.cpp cfiles/mainrand.h
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

cfiles/mainrand.o: cfiles/mainrand.cpp cfiles/mainrand.h
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
