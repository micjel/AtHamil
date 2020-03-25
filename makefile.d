obj/AtLibrary.o : src/AtLibrary.f90 
obj/AtomicHamiltonian.o : src/AtomicHamiltonian.F90 obj/TwoBodySpinDipole.o obj/TwoBodySpinOrbit.o obj/TwoBodyOrbitOrbit.o obj/TwoBodySpinContact.o obj/TwoBodyDarwin.o obj/TwoBodyCoulomb.o obj/AtLibrary.o obj/OneBodyTerms.o obj/EleSingleParticleState.o obj/ElectronTwoBodySpace.o obj/LinAlgLib.o 
obj/EleSingleParticleState.o : src/EleSingleParticleState.F90 obj/AtLibrary.o 
obj/ElectronTwoBodySpace.o : src/ElectronTwoBodySpace.F90 obj/AtLibrary.o obj/EleSingleParticleState.o 
obj/OneBodyTerms.o : src/OneBodyTerms.F90 obj/AtLibrary.o obj/EleSingleParticleState.o obj/LinAlgLib.o 
obj/TwoBodyCoulomb.o : src/TwoBodyCoulomb.F90 obj/AtLibrary.o obj/ElectronTwoBodySpace.o obj/OneBodyTerms.o 
obj/TwoBodyDarwin.o : src/TwoBodyDarwin.F90 obj/AtLibrary.o obj/ElectronTwoBodySpace.o obj/OneBodyTerms.o 
obj/TwoBodyOrbitOrbit.o : src/TwoBodyOrbitOrbit.F90 obj/AtLibrary.o obj/ElectronTwoBodySpace.o obj/OneBodyTerms.o 
obj/TwoBodySpinContact.o : src/TwoBodySpinContact.F90 obj/AtLibrary.o obj/ElectronTwoBodySpace.o obj/OneBodyTerms.o 
obj/TwoBodySpinDipole.o : src/TwoBodySpinDipole.F90 obj/AtLibrary.o obj/ElectronTwoBodySpace.o obj/OneBodyTerms.o 
obj/TwoBodySpinOrbit.o : src/TwoBodySpinOrbit.F90 obj/AtLibrary.o obj/ElectronTwoBodySpace.o obj/OneBodyTerms.o 
obj/ClassSys.o : main/ClassSys.f90 
obj/AtHamilInput.o : main/AtHamilInput.F90 
obj/AtHamilMain.o : main/AtHamilMain.F90 obj/AtHamilInput.o obj/AtomicHamiltonian.o obj/EleSingleParticleState.o obj/AtLibrary.o 
obj/LinAlgLib.o : submodules/LinAlgf90/src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatVecSingle.o obj/MatrixComplex.o obj/MatrixDouble.o obj/MatrixSingle.o obj/VectorComplex.o obj/VectorDouble.o obj/VectorSingle.o obj/SingleDoubleComplex.o obj/LinAlgParameters.o 
obj/LinAlgParameters.o : submodules/LinAlgf90/src/LinAlgParameters.f90 
obj/MatVecComplex.o : submodules/LinAlgf90/src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MatVecDouble.o : submodules/LinAlgf90/src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatVecSingle.o : submodules/LinAlgf90/src/MatVecSingle.f90 obj/MatrixSingle.o obj/VectorSingle.o obj/LinAlgParameters.o 
obj/MatrixComplex.o : submodules/LinAlgf90/src/MatrixComplex.f90 obj/VectorComplex.o obj/LinAlgParameters.o 
obj/MatrixDouble.o : submodules/LinAlgf90/src/MatrixDouble.f90 obj/VectorDouble.o obj/LinAlgParameters.o 
obj/MatrixSingle.o : submodules/LinAlgf90/src/MatrixSingle.f90 obj/VectorSingle.o obj/LinAlgParameters.o 
obj/SingleDoubleComplex.o : submodules/LinAlgf90/src/SingleDoubleComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/MatrixDouble.o obj/VectorDouble.o obj/MatrixSingle.o obj/VectorSingle.o 
obj/VectorComplex.o : submodules/LinAlgf90/src/VectorComplex.f90 obj/LinAlgParameters.o 
obj/VectorDouble.o : submodules/LinAlgf90/src/VectorDouble.f90 obj/LinAlgParameters.o 
obj/VectorSingle.o : submodules/LinAlgf90/src/VectorSingle.f90 obj/LinAlgParameters.o 
