
INC_DIRS = /opt/homebrew/opt/libomp/include
LIB_DIRS = /opt/homebrew/opt/libomp/lib
LIBS = -lomp

CXX = c++
CXXFLAGS = -Xpreprocessor -fopenmp -march=native -std=c++17

HDRS := crs.hpp precond.hpp
SRCS := crs_io.cpp cg.cpp

TARGET = cg dcg iccg

all: $(TARGET)

clean:
	rm -f $(TARGET) *.o

cg: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -I$(INC_DIRS) $(SRCS) -L$(LIB_DIRS) $(LIBS) -o $@

dcg: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -DJAC -I$(INC_DIRS) $(SRCS) -L$(LIB_DIRS) $(LIBS) -o $@

iccg: $(SRCS) ic0.cpp $(HDRS) ic0.hpp
	$(CXX) $(CXXFLAGS) -DICP -I$(INC_DIRS) $(SRCS) ic0.cpp -L$(LIB_DIRS) $(LIBS) -o $@

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<