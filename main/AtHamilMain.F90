program AtHamil
  use AtLibrary, only: isfile
  use EleSingleParticleState
  use AtomicHamiltonian
  use AtomicOperator
  use AtHamilInput, only: InputParameters
  type(InputParameters) :: params
  type(EleOrbits) :: sps
  type(EleTwoBodySpace) :: ms
  type(AtomicHamil) :: h
  type(Op) :: ops
  character(:), allocatable :: f
  call params%init()

  call sps%init(params%emax, lmax=params%lmax, mode=params%basis)
  call ms%init(sps, params%e2max, params%zeta)

  select case(params%OpName)
  case("Coulomb", "Breit", "Breit_1b")
    call h%init(ms, params%OpName)
    f = trim(params%file_name)
    if(f == "default") f = h%GetFileName()
    if(isfile(f)) then
      write(*,"(a)") trim(f), " already exsits."
      stop
    end if
    call h%set(params%NMesh, params%rmax, params%NMesh_Mom, params%pmax)
    call h%writef(f)
    call h%fin()
  case("MassShift","NormalMassShift","SpecificMassShift","FieldShift")
    call ops%init(params%OpName, ms, 0, 1)
    call ops%SetOp(params%NMesh, params%pmax)
    call ops%WriteOp(params%file_name)
    call ops%fin()
  case default

  end select
  call ms%fin()
  call sps%fin()
end program AtHamil
