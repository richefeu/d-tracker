# At least g77 OR gfortran must be installed!
CC = $(shell if [ `which g77` ]; then echo "g77"; else  echo "gfortran"; fi)
#CC = gfortran
#CC = g77

FLAGS = -fno-backslash

# Add here the list of source files
SRC = BilinearInterpol.f pearsn.f string.f

# Do not take care about the following
OBJS = $(subst .f,.o,$(SRC))

.SUFFIXES: .f
.f.o: $(SRC)
	$(CC)  $(FLAGS) -c $<

all: $(OBJS)

clean:
	rm -f *.o 

