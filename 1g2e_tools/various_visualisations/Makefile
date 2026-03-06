APPS = LocalisationDeformation DessinEch Deplacement-2 RotationMoyenne draw-real8-t
LIB_DIR = ../Library

# At least g77 OR gfortran must be installed!
CC = $(shell if [ `which g77` ]; then echo "g77"; else echo "gfortran"; fi)

FLAGS = -fno-backslash

all: $(APPS)

LocalisationDeformation: LocalisationDeformation.f Makefile
	$(CC) $(FLAGS) -c $@.f
	$(CC) $(FLAG) -o $@.e $@.o -L/usr/local/lib -lpgplot

DessinEch: DessinEch.f Makefile
	$(CC) $(FLAGS) -c $@.f
	$(CC) $(FLAGS) -o $@.e $@.o -L/usr/local/lib -lpgplot

Deplacement-2: Deplacement-2.f Makefile
	$(CC) $(FLAGS) -c $@.f
	$(CC) $(FLAGS) -o $@.e $@.o $(LIB_DIR)/BilinearInterpol.o -L/usr/local/lib -lpgplot

RotationMoyenne: RotationMoyenne.f Makefile
	$(CC) $(FLAGS) -c $@.f
	$(CC) $(FLAGS) -o $@.e $@.o $(LIB_DIR)/string.o 

draw-real8-t: draw-real8-t.f Makefile
	$(CC) $(FLAGS) -c $@.f -fbackslash -ffixed-line-length-132
	$(CC) $(FLAGS) -o $@.e $@.o -L/usr/local/lib  -lpgplot 

clean:
	rm -f *.e *.o 

