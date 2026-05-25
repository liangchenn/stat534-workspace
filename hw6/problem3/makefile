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

regmodels.o: regmodels.cpp matrices.h
	gcc -g -c regmodels.cpp -o regmodels.o ${MATRICES_INCLUDE}

matrices: main.o matrices.o regmodels.o
	gcc main.o matrices.o regmodels.o -o matrices ${MATRICES_LIB} ${MATRICES_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MATRICES_OBJS}
