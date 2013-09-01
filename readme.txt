This program is currently designed to be compiled for quantum espresso 4.2.1

Instructions to get up and running:
1. Download Quantum Espresso program and examples and decompress and change the current directory to the Quantum Espresso directory
2. configure the build to disable MPI and enable OpenMP ./configure --disable-parallel --enable-openmp
3. Make a folder called kubocalc and make it the current directory
4. copy kubocalc.f90 and Makefile into folder
5. type make - this will build the kubocalc.x executable
6. if you want to run MPI executables, reconfigure and rebuild:
   cd ..
   make clean
   ./configure
   make all
7. copy the kubocalc.x executable into the bin directory
   cp kubocalc/kubocalc.x BIN/

You now have a working installation and are ready to run kubo-greenwood calculations

There is an example for solid aluminum:
To run the 4 atom Aluminum unit cell calculation, make a directory in you you need the following files
run_md_Al4                 - this file runs the molecular dynamics calculation to get an atomic configuration
run_scf_Al4                - this file runs the SCF calculation for that atomic configuration with higher precision
run_kubo_Al4               - this file runs the kubo-greenwood calculation. The output file is in results\optcond.dat

of course, this gives only a single data point and you have to do an ensemble calculation.....
to make things a little easier, you can run
run_all_Al4                - this file does all the work. It produces an number of intermediate files you can use to check your work, but the final output is
                           - optcondAverage_nat${nat}_k${k}.txt where ${nat} is the number of atoms (in this case 4) and ${k} is the number of k points for the SCF calculation

These files should be very easy to adjust to your particular case. As examples of how to do this there are other file sets. First,
run_md_Al32                - this file runs the molecular dynamics calculation to get an atomic configuration
run_scf_Al32               - this file runs the SCF calculation for that atomic configuration with higher precision
run_kubo_Al32              - this file runs the kubo-greenwood calculation. The output file is in results\optcond.dat
run_all_Al32               - this file does all the work.

And also ongoing work for Magnesium Dicilicide:
run_md_Mg2Si24
run_scf_Mg2Si24
run_kubo_Mg2Si24
run_all_Mg2Si
