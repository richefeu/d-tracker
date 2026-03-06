APPS = TemporalFlucDispCorrelation MakeFluctuationMovie AutoCorrelSpace-couple StatVelocities
LIB_DIR = ../Library

# At least g77 OR gfortran must be installed!
CC = $(shell if [ `which g77` ]; then echo "g77"; else echo "gfortran"; fi)

all: $(APPS)

TemporalFlucDispCorrelation: TemporalFlucDispCorrelation.f Makefile
	$(CC) -c $@.f -O3
	$(CC) -o $@.e $@.o $(LIB_DIR)/BilinearInterpol.o $(LIB_DIR)/pearsn.o 

MakeFluctuationMovie: MakeFluctuationMovie.f Makefile
	$(CC) -c $@.f -O3
	$(CC) -o $@.e $@.o $(LIB_DIR)/BilinearInterpol.o -lpgplot

AutoCorrelSpace-couple: AutoCorrelSpace-couple.f Makefile
	$(CC) -c $@.f -O3
	$(CC) -o $@.e $@.o $(LIB_DIR)/BilinearInterpol.o 


StatVelocities: StatVelocities.f Makefile
	$(CC) -c $@.f -O3
	$(CC) -o $@.e $@.o $(LIB_DIR)/BilinearInterpol.o  


clean:
	rm -f *.e *.o 

