# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -O2 -std=c++20 -fopenmp

# Debug flag
ifeq ($(DEBUG),1)
    CXXFLAGS += -g
endif

# Source files
CPP_SRCS = main.cpp parsim.cpp initParticles.cpp parsimUtils.cpp

# Header files
HEADERS = ds.hpp initParticles.hpp parsim.hpp parsimUtils.hpp pointCalculations.hpp

# Object files
CPP_OBJS = $(CPP_SRCS:.cpp=.o)

# Executable name
TARGET = parsim

# Default rule
all: $(TARGET)

# Build executable
$(TARGET): $(CPP_OBJS)
	$(CXX) $(CXXFLAGS) $(CPP_OBJS) -o $(TARGET)

# Compile C++ source files into object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f $(CPP_OBJS) $(TARGET)

# Phony targets
.PHONY: all clean
