CC = g++
LDIR = ../../lib
IDIR = ../../include
CFLAGS = -I$(IDIR) -L$(LDIR)
LIBS = -lgsl -lgslcblas
OBJ = mainrand.o cycledef.o main.o

gd_random.out: $(OBJ)
	$(CC) -o $@ $^ -shared $(CFLAGS) $(LIBS)

main.o: main.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

cycledef.o: cycledef.cpp mainrand.h
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

mainrand.o: mainrand.cpp mainrand.h
	$(CC) -c -o $@ $< $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
