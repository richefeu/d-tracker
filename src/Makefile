APP = tracker

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# You can install clang-omp by using homebrew: brew install clang-omp
	CC = g++-15
	CC_OPT = -fopenmp -O3 -W -Wall -std=c++17
else
	CC = g++
	CC_OPT = -fopenmp -O3 -pedantic -W -Wall -Wno-deprecated -std=c++0x
endif

# Ensure the toofus directory exists
ifeq ($(wildcard ./toofus),)
$(info Cloning ToOfUs repository)
$(shell git clone https://github.com/richefeu/toofus.git > /dev/null 2>&1)
endif

# You have to compile libRaw (available in the folder ThirdParty) and install it
# ./configure, make, sudo make install
# You can also install it with Homebrew (on apple-mac), or something else (apt-get, etc.).
CONFIG_LIBRAW = `pkg-config --cflags --libs libraw`
CONFIG_LIBTIFF = `pkg-config --cflags --libs libtiff-4`

.PHONY: all clean clone_toofus

# All the files that are part of the source code
SOURCES = $(APP).cpp $(APP).hpp thumbnail.hpp 

all: $(APP)

$(APP): run.cpp $(SOURCES) Makefile
	$(CC) $(CC_OPT) -o $(APP) run.cpp $(CONFIG_LIBRAW) $(CONFIG_LIBTIFF) 

clean:
	rm -f *.o $(APP)  

clone_toofus:
	@if [ ! -d "./toofus" ]; then \
		@echo "\033[0;32m-> CLONING ToOfUs (Tools often used)\033[0m"; \
		git clone https://github.com/richefeu/toofus.git; \
	fi