

CC=g++
OBJ= selfe.o
LIBDIR=/usr/local/lib
CFLAGS=-std=c++14 -I.:/usr/local/include -L$(LIBDIR) -lgsl -lgslcblas -larmadillo
DEPS=  selfe.h structs.h declarations.h
selfe: $(OBJ) 
	$(CC)  -o $@  $(OBJ) $(CFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o
