EBROOTGSL ?= $(shell pkg-config --variable=prefix gsl 2>/dev/null)
MATRICES_LSTDFLG = -lstdc++ -lgsl -lgslcblas -lm
MATRICES_INCLUDE = -I${EBROOTGSL}/include
MATRICES_LIB = -L${EBROOTGSL}/lib
MATRICES_OBJS = matrices

all: ${MATRICES_OBJS}

matrices.o: matrices.cpp matrices.h
	gcc -g -c matrices.cpp -o matrices.o ${MATRICES_INCLUDE}

main.o: main.cpp matrices.h
	gcc -g -c main.cpp -o main.o ${MATRICES_INCLUDE}

matrices: main.o matrices.o
	gcc main.o matrices.o -o matrices ${MATRICES_LIB} ${MATRICES_LSTDFLG} && ./matrices

clean:
	rm -f *.o
	rm -f ${MATRICES_OBJS}
