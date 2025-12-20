
INC_DIRS = /opt/homebrew/opt/libomp/include
LIB_DIRS = /opt/homebrew/opt/libomp/lib
LIBS = -lomp

CXX = c++
CXXFLAGS = -Xpreprocessor -fopenmp -I$(INC_DIRS)

HDRS := crs.hpp
SRCS := crs_io.cpp cg.cpp
OBJS := $(SRCS:.cpp=.o)

TARGET = cg

all: $(TARGET)

clean:
	rm -f $(TARGET) *.o

$(TARGET): $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(OBJS) -L$(LIB_DIRS) $(LIBS) -o $@

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<