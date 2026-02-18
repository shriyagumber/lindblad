SRC := src/solve_me2.f95 \
       src/solve_tdse.f95 \
       src/readHnL.f95 \
       src/printfssh.f95 \
       src/inter.f95 \
       src/fssh.f95 \
       src/main_lindblad.f95

OBJ := $(patsubst src/%.f95,%.o,$(SRC))

all: run_lindblad

run_lindblad: $(OBJ)
	gfortran $(OBJ) -o $@

%.o: src/%.f95
	gfortran -c $<

clean:
	rm -f *.o *.mod run_lindblad
