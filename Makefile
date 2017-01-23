CXX = g++
CXXFLAGS = -Wall -fexceptions -O2 -std=c++11 

PRG = introgress
CLASSES = Indiv.o Maths.o Network.o Population.o
OBJ = main.o $(CLASSES)

all: $(PRG)

$(PRG): $(OBJ)
	$(CXX) $(CXXFLAGS) -larmadillo -o $@ $(OBJ)

.cpp.o:
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(PRG) $(OBJ) $(TEST)
