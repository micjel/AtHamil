# AtHamil
This code generates the matrix elements of atomic Hamiltonians with the Hamonic-Oscillator, Hydrogen, Laguerre-function basis sets.
For the details of Laguerre-function basis, see A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).

Currently, it can generate only non-relativistic Coulomb Hamiltonian.

Python code is only for the quick test.

## Requirements
- Fortran compiler
- BLAS and LAPACK installation
- GNU scientific library
- zlib

On Cedar, it is suggested to install the modules using 
```
module load gcc gsl flexiblas
```
and the default configuration in the in the Makefile corresponds to this setup. 

## Installation
1. Download the code with the LinAlgf90 submodule using
    ```
    git clone --recursive https://github.com/micjel/AtHamil.git
    ``` 
    into the directory you wish to store AtHamil in. Configure the Makefile in the main directory to work with your installation if you are not running on Oak or Cedar.  

2. Once the build parameters are configured correctly, run `make` in the terminal. If the installation is successful, the file AtHamil.exe should be updated in exe/.

3. Run `make install` to copy the .exe to ~/bin; ensure that this is added to your path. Now the job script AtomicHamil.py can be used to submit jobs to create the Hamiltonians or run locally.

## Known issues
- Compiling AtHamil with OpenBLAS instead of FlexiBLAS on Cedar leads to AtHamil.exe hanging during the execution. I am unsure why this occurs.
- If you are compiling with gcc 10 and above, you should use `DEBUG_MODE=off` in line 18; this causes the type mismatches in the submodule to be correctly treated as warnings.
- If you are compiling with gcc version 9 and below, delete the flag `-fallow-argument-mismatch` from the compilation.