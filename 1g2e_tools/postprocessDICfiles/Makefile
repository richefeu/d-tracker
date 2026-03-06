APP = PostProcessDIC

# At least g77 OR gfortran must be installed!
CC = $(shell if [ `which g77` ]; then echo "g77"; else echo "gfortran"; fi)


SOURCES = $(APP).f

all: $(APP)

$(APP) : $(SOURCES) Makefile
	$(CC) -o $(APP).e $(APP).f -fno-backslash -O3

clean:
	rm -f $(APP).e 

