APP = tracker
GUI = tracker-gui
LIB = libtracker.a

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# You can install clang-omp by using homebrew: brew install clang-omp
	CC = g++-15
	CC_OPT = -fopenmp -O3 -W -Wall -Wextra -Wpedantic -std=c++17
else
	CC = g++
	CC_OPT = -fopenmp -O3 -pedantic -W -Wall -Wno-deprecated -std=c++17
endif

# Ensure the toofus directory exists
ifeq ($(wildcard ./toofus),)
$(info Cloning ToOfUs repository)
$(shell git clone https://github.com/richefeu/toofus.git > /dev/null 2>&1)
endif

# You have to compile libRaw (available in the folder ThirdParty) and install it
# ./configure, make, sudo make install
# You can also install it with Homebrew (on apple-mac), or something else (apt-get, etc.).
CONFIG_LIBRAW = $(shell pkg-config --cflags --libs libraw | sed 's/-lstdc++//g')
CONFIG_LIBTIFF = $(shell pkg-config --cflags --libs libtiff-4)

# raylib (for the GUI application tracker-gui). Installed e.g. with: brew install raylib
CONFIG_RAYLIB = $(shell pkg-config --cflags --libs raylib)
ifeq ($(UNAME_S),Darwin)
	# Frameworks required by raylib on macOS
	CONFIG_RAYLIB += -framework Cocoa -framework IOKit -framework OpenGL
endif

.PHONY: all clean clone_toofus get_raygui

# All the files that are part of the core source code (compiled into libtracker.a)
SOURCES = $(APP).cpp $(APP).hpp thumbnail.hpp interpol.hpp

all: $(APP) $(GUI)

# --- Static library that holds the whole tracker core -------------------------
$(LIB): $(SOURCES) Makefile
	$(CC) $(CC_OPT) -c $(APP).cpp -o $(APP).o $(CONFIG_LIBRAW) $(CONFIG_LIBTIFF)
	ar rcs $(LIB) $(APP).o

# --- Command-line front-end ---------------------------------------------------
$(APP): run.cpp $(LIB) Makefile
	$(CC) $(CC_OPT) -o $(APP) run.cpp $(LIB) $(CONFIG_LIBRAW) $(CONFIG_LIBTIFF)

# --- raylib/raygui GUI front-end ----------------------------------------------
# raygui.h (single-header) declenche des -Wmissing-field-initializers : on les
# neutralise uniquement pour cette cible.
$(GUI): tracker_gui.cpp $(LIB) raygui.h guifont_ttf.h Makefile
	$(CC) $(CC_OPT) -Wno-missing-field-initializers -Wno-enum-compare -Wno-unused-parameter -o $(GUI) tracker_gui.cpp $(LIB) $(CONFIG_RAYLIB) $(CONFIG_LIBRAW) $(CONFIG_LIBTIFF)

# raygui is a single-header UI library for raylib. Fetch it if missing.
raygui.h:
	@echo "\033[0;32m-> Downloading raygui.h\033[0m"
	curl -fsSL -o raygui.h https://raw.githubusercontent.com/raysan5/raygui/4.0/src/raygui.h

# Fonte d'interface (Open Sans, licence OFL) embarquee dans le binaire.
# On telecharge le TTF puis on le transforme en tableau C avec 'xxd -i'.
guifont.ttf:
	@echo "\033[0;32m-> Downloading UI font (Open Sans)\033[0m"
	curl -fsSL -o guifont.ttf https://github.com/googlefonts/opensans/raw/main/fonts/ttf/OpenSans-Regular.ttf

guifont_ttf.h: guifont.ttf
	@echo "\033[0;32m-> Embedding UI font into guifont_ttf.h\033[0m"
	xxd -i guifont.ttf > guifont_ttf.h

get_raygui: raygui.h

clean:
	rm -f *.o $(APP) $(GUI) $(LIB)

clone_toofus:
	@if [ ! -d "./toofus" ]; then \
		@echo "\033[0;32m-> CLONING ToOfUs (Tools often used)\033[0m"; \
		git clone https://github.com/richefeu/toofus.git; \
	fi