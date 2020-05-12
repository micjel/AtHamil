module AtHamilInput
  implicit none

  public :: InputParameters
  private :: InitInputParameters
  private :: CopyInputParameters
  private :: PrintInputParameters

  type :: InputParameters
    ! Input parameters
    ! General
    real(8) :: zeta ! length scale of basis funcgtion, in unit of (a.u.)^-1
    character(:), allocatable :: basis ! Basis function HO, AO, LO
    integer :: NMesh
    real(8) :: rmax
    integer :: NMesh_Mom
    real(8) :: pmax
    character(:), allocatable :: file_name
    character(:), allocatable :: opname
    integer :: emax, e2max, lmax
    logical :: test_mode
    logical :: count_memory
  contains
    procedure :: init => InitInputParameters
    procedure :: CopyInputParameters
    generic :: assignment(=) => CopyInputParameters
    procedure :: prt => PrintInputParameters
  end type InputParameters
contains
  subroutine InitInputParameters(params)
    class(InputParameters), intent(inout) :: params
    real(8) :: zeta = 1.d0
    integer :: NMesh = 400
    real(8) :: rmax=1.d3
    integer :: NMesh_Mom = 400
    real(8) :: pmax=1.d3
    character(256) :: basis = "LO"
    character(512) :: file_name = "default"
    character(256) :: opname = "Breit"
    integer :: emax = 6
    integer :: e2max=12
    integer :: lmax = -1
    logical :: test_mode = .false.
    logical :: count_memory = .false.
    logical :: ex

    ! input file name
    character(512) :: inputfile
    namelist /input/ zeta, basis, file_name, emax, e2max, lmax, test_mode, count_memory, &
        & rmax, NMesh, pmax, NMesh_Mom, opname
    call getarg(1, inputfile)
    inquire(file = inputfile, exist = ex)
    if(.not.ex) then
      write(*,'(3a)') 'Error: Input file: ', trim(inputfile), ' was not found.'
      stop
    end if
    open(118, file = inputfile, status = 'old')
    read(118, nml=input)
    close(118)
    params%file_name = file_name
    params%emax = emax
    params%lmax = lmax
    params%e2max = e2max
    params%test_mode = test_mode
    params%zeta = zeta
    params%basis = basis
    params%rmax = rmax
    params%Nmesh = NMesh
    params%pmax = pmax
    params%Nmesh_Mom = NMesh_Mom
    params%opname = opname
    if(params%lmax==-1) params%lmax = params%emax

    params%count_memory = count_memory
  end subroutine InitInputParameters

  subroutine CopyInputParameters(params2, params1)
    class(InputParameters), intent(inout) :: params2
    type(InputParameters), intent(in) :: params1


    params2%file_name= params1%file_name
    params2%emax     = params1%emax
    params2%lmax     = params1%lmax
    params2%e2max    = params1%e2max
    params2%test_mode= params1%test_mode
    params2%zeta     = params1%zeta
    params2%basis    = params1%basis
    params2%opname   = params1%opname
  end subroutine CopyInputParameters

  subroutine PrintInputParameters(params,unt)
    class(InputParameters), intent(in) :: params
    integer, intent(in), optional :: unt
    integer :: iunit
    if(present(unt)) then
      iunit = unt
    else
      iunit = 6
    end if
    write(iunit,*)
    write(iunit,'(a)') '################### Input Parameters #####################'
#ifdef VERSION
    write(*,"(2a)") "# NuHamil version is ", trim(VERSION)
#endif
    write(iunit,"(a)") "# Calculation for the matrix elements of Atomic Hamiltonian"
    write(iunit,"(2a)") "# Operator", trim(params%opname)
    write(iunit,'(a, i3, a, i3, a, i3)') '# emax: ', params%emax, ', e2max: ', params%e2max, &
        & ", lmax: ", params%lmax
    write(iunit,'(3a,f6.2)') '# basis: ', trim(params%basis), ", basis parameter: ", params%zeta
    write(iunit,*)
  end subroutine PrintInputParameters
end module AtHamilInput
