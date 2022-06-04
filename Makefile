CC = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp
BOOSTPATH = ../cpp/boost_1_75_0
INCLUDE = -I $(BOOSTPATH) -I ./include/
LIB = -L $(BOOSTPATH)
SRC = ./src/
TEST = 

CPPFILES = $(SRC)main.cpp \
				$(SRC)ENO3_1d.cpp \
				$(SRC)ENO3_2d.cpp \
				$(SRC)SemiLagrangian_Cancer.cpp \

main: $(CPPFILES)
	$(CC) $(CXXFLAGS) -o $@ $^ $(INCLUDE) $(LIB)

run: main
		./main $(TEST)

debug: $(CPPFILES)
	g++ -o1 -g -std=c++17 $(INCLUDE) $^ $(LIB)
	gdb -ex=r --args 1 $(TEST)

clean:
	rm *.dat main
