EBROOTGSL ?= $(shell pkg-config --variable=prefix gsl 2>/dev/null)
CC = gcc
MATRICES_LSTDFLG = -lstdc++ -lgsl -lgslcblas -lm
MATRICES_INCLUDE = -I${EBROOTGSL}/include
MATRICES_LIB = -L${EBROOTGSL}/lib

SRC = main.cpp matrices.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = matrices
EXAMPLE_SRC = main-example.cpp
EXAMPLE_OBJ = $(EXAMPLE_SRC:.cpp=.o)
EXAMPLE_TARGET = main-example
CHECK_SRC = check.cpp
CHECK_OBJ = $(CHECK_SRC:.cpp=.o)
CHECK_TARGET = check-cov

.PHONY: all clean run test check

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(MATRICES_LIB) $(MATRICES_LSTDFLG)

$(EXAMPLE_TARGET): $(EXAMPLE_OBJ)
	$(CC) $(EXAMPLE_OBJ) -o $(EXAMPLE_TARGET) $(MATRICES_LIB) $(MATRICES_LSTDFLG)

$(CHECK_TARGET): $(CHECK_OBJ)
	$(CC) $(CHECK_OBJ) -o $(CHECK_TARGET) $(MATRICES_LIB) $(MATRICES_LSTDFLG)


%.o: %.cpp matrices.h
	$(CC) $(CFLAGS) $(MATRICES_INCLUDE) -c $< -o $@


run: $(TARGET)
	./$(TARGET)

test: $(EXAMPLE_TARGET)
	./$(EXAMPLE_TARGET)

check: $(CHECK_TARGET)
	./$(CHECK_TARGET)

clean:
	rm -f $(OBJ) $(TARGET) $(EXAMPLE_OBJ) $(EXAMPLE_TARGET) $(CHECK_OBJ) $(CHECK_TARGET)

print:
	echo $(CC)
