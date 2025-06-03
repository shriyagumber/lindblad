gfortran -c solve_me2.f95
gfortran -c solve_tdse.f95
gfortran -c readHnL.f95
gfortran -c printfssh.f95
gfortran -c inter.f95
gfortran -c fssh.f95
gfortran -c main_lindblad.f95

gfortran readHnL.o fssh.o inter.o solve_me2.o solve_tdse.o printfssh.o main_lindblad.o -o run_lindblad

rm *.o *.mod
