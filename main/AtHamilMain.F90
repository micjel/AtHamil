program AtHamil
  use AtLibrary, only: isfile
  use EleSingleParticleState
  use AtomicHamiltonian
  use AtHamilInput, only: InputParameters
  type(InputParameters) :: params
  type(EleOrbits) :: sps
  type(EleTwoBodySpace) :: ms
  type(AtomicHamil) :: h
  character(:), allocatable :: f
  call params%init()

  call sps%init(params%emax, lmax=params%lmax, mode=params%basis)
  call ms%init(sps, params%e2max, params%zeta)

  call h%init(ms)
  f = trim(params%file_name)
  if(f == "default") f = h%GetFileName()
  if(isfile(f)) then
    write(*,"(a)") trim(f), " already exsits."
    return
  end if
  call h%set(params%NMesh, params%rmax)
  call h%writef(f)

  call h%fin()
  call ms%fin()
  call sps%fin()
end program AtHamil
