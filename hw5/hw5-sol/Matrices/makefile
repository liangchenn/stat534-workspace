EBROOTOPENBLAS ?= /opt/homebrew/opt/openblas
MATRICES_LSTDFLG = -lstdc++ -lopenblas -lm
MATRICES_INCLUDE = -I${EBROOTOPENBLAS}/include
MATRICES_LIB = -L${EBROOTOPENBLAS}/lib
MATRICES_OBJS = matrices

all: ${MATRICES_OBJS}

matrices.o: matrices.cpp matrices.h
	gcc -g -c matrices.cpp -o matrices.o ${MATRICES_INCLUDE}

main.o: main.cpp matrices.h
	gcc -g -c main.cpp -o main.o ${MATRICES_INCLUDE}

matrices: main.o matrices.o
	gcc main.o matrices.o -o matrices ${MATRICES_LIB} ${MATRICES_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MATRICES_OBJS}
