
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

all: seeCONF # extractBoxData

clean:
	rm -f *.o seeCONF

# The application that visualize the conf files
seeCONF: seeCONF.cpp seeCONF.hpp confProcess.hpp confProcess.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -Wno-deprecated $(GLLINK)

extractBoxData: extractBoxData.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
