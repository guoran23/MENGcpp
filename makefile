# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

# Target executable
TARGET = simulation

# Source files
SRCS = main.cpp Field.cpp Particle.cpp Equilibrium.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Build rules
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

