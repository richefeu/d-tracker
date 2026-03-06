APP = tracker
#SVNDEV := -D'SVN_REV="$(shell svnversion -n .)"'
SVNDEV := -D'SVN_REV="SVN_Rev_Deprecated"'

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# You can install clang-omp by using homebrew: brew install clang-omp
	CC = g++-14
	CC_OPT = -fopenmp -O3 -W -Wall -std=c++11 $(SVNDEV)  
else
	CC = g++
	CC_OPT = -fopenmp -O3 -pedantic -W -Wall -Wno-deprecated -std=c++0x $(SVNDEV)
endif


# You have to compile libRaw (available in the folder ThirdParty) and install it
# ./configure, make, sudo make install
# You can also install it with Macports or Homebrew (on apple-mac), or something else (apt-get, etc.).
CONFIG_LIBRAW = `pkg-config --cflags --libs libraw`
CONFIG_LIBTIFF = `pkg-config --cflags --libs libtiff-4`


# The following overwrites all previous (automatic) settings
ifneq ($(wildcard ./local_config.inc),) 
	include ./local_config.inc
endif
# The file local_config.inc may contain something like this:
#CC = g++-4.9
#CC_OPT = -fopenmp -O3 -pedantic -W -Wall -Wextra -Wno-long-long $(SVNDEV) -std=c++0x -I/usr/local/include


# All the files that are part of the source code
SOURCES = $(APP).cpp $(APP).hpp \
thumbnail.hpp \
ColorTable.hpp ColorTable.cpp \
powell_nr3.hpp \
distortion_correction.cpp \
pattern_quality.cpp \
subpixel_centers.cpp \
synthetic_images.cpp \
post_process.cpp \
visu_process.cpp \
grayLevelAnalysis.cpp

all: $(APP)

$(APP): $(SOURCES) $(APP)SH.cpp Makefile
	$(CC) $(CC_OPT) -o $(APP) $(APP)SH.cpp $(CONFIG_LIBRAW) $(CONFIG_LIBTIFF) 

mapNeighbor:mapNeighbor.cpp
	$(CC) $(CC_OPT) -o mapNeighbor mapNeighbor.cpp $(CONFIG_LIBTIFF) $(CONFIG_LIBRAW) 
	
lineCorrel:lineCorrel.cpp
	$(CC) $(CC_OPT) -o lineCorrel lineCorrel.cpp $(CONFIG_LIBTIFF) $(CONFIG_LIBRAW) 

clean:
	rm -f *.o $(APP) mapNeighbor lineCorrel

depend:
	makedepend -Y $(SOURCES)

# DO NOT DELETE

tracker.o: tracker.hpp ColorTable.cpp ColorTable.hpp powell_nr3.hpp
tracker.o: interpol.cpp delaunay2D.hpp thumbnail.hpp image_io.cpp
tracker.o: distortion_correction.cpp pattern_quality.cpp post_process.cpp
tracker.o: subpixel_centers.cpp synthetic_images.cpp visu_process.cpp
tracker.o: ColorTable.cpp ColorTable.hpp powell_nr3.hpp interpol.cpp
tracker.o: delaunay2D.hpp thumbnail.hpp image_io.cpp
ColorTable.o: ColorTable.hpp
