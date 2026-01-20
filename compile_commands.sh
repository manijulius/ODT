!#/bin/bash


# edit Labinit.f file to configure the simulation initial and boundary conditions
# compile and execute Labinit to produce input argument file "LabExppar.dat"
gfortran -O3 Labinit.f  -o Labinit.exe

# execute Labinit to produce input argument file "LabExppar.dat" for LabExp
./Labinit.exe










#compile main ODT program
gfortran -O3 LabExp.f  -o LabExp.exe

# Run the executable with experiment_name "test" and random perturbation seed 001
# simulation output will be stored under output/test
nohup ./LabExp.exe steady_state test 001 &> 001.log &


