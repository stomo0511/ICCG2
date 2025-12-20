
INC_DIRS = /opt/homebrew/opt/libomp/include
LIB_DIRS = /opt/homebrew/opt/libomp/lib
LIBS = -lomp

CXX = c++
CXXFLAGS = -Xpreprocessor -fopenmp -I$(INC_DIRS)

HDRS := crs.hpp
SRCS := crs_io.cpp cg.cpp

TARGET = cg dcg

all: $(TARGET)

clean:
	rm -f $(TARGET) *.o

cg: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(SRCS) -L$(LIB_DIRS) $(LIBS) -o $@

dcg: $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -DJAC $(SRCS) -L$(LIB_DIRS) $(LIBS) -o $@

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<