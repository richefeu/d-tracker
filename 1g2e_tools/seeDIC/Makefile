
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  # The compiler to be used
  CXX = g++-15
  # The list of flags passed to the compiler
  CXXFLAGS = -O3 -Wall -std=c++11
  # Link flags for OpenGL and glut
  GLLINK = `pkg-config --cflags --libs gl glu glut`
else
  # The compiler to be used
  CXX = g++
  # The list of flags passed to the compiler
  CXXFLAGS = -O3 -Wall -std=c++11
  # Link flags for OpenGL and glut
  GLLINK = -lGLU -lGL -L/usr/X11R6/lib -lglut -lXmu -lXext -lX11 -lXi
endif

.PHONY: all clean

all: seeDIC

clean:
	rm -f *.o seeDIC

# The application that visualize the dic files
seeDIC: seeDIC.cpp seeDIC.hpp dicProcess.hpp dicProcess.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -Wno-deprecated $(GLLINK)


