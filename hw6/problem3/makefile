EBROOTGSL ?= $(shell pkg-config --variable=prefix gsl 2>/dev/null)
REGRESSION_LSTDFLG = -lstdc++ -lgsl -lgslcblas -lm
REGRESSION_INCLUDE = -I${EBROOTGSL}/include
REGRESSION_LIB = -L${EBROOTGSL}/lib
REGRESSION_OBJS = regression

all: ${REGRESSION_OBJS}

matrices.o: matrices.cpp matrices.h
	gcc -g -c matrices.cpp -o matrices.o ${REGRESSION_INCLUDE}

main.o: main.cpp matrices.h
	gcc -g -c main.cpp -o main.o ${REGRESSION_INCLUDE}

regmodels.o: regmodels.cpp matrices.h
	gcc -g -c regmodels.cpp -o regmodels.o ${REGRESSION_INCLUDE}

regression: main.o matrices.o regmodels.o
	gcc main.o matrices.o regmodels.o -o regression ${REGRESSION_LIB} ${REGRESSION_LSTDFLG} && ./regression

clean:
	rm -f *.o
	rm -f ${REGRESSION_OBJS}
