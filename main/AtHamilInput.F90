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
    character(:), allocatable :: file_name
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
    real(8) :: zeta = 4.d0
    character(256) :: basis = "LO"
    character(512) :: file_name = "default"
    integer :: emax = 6
    integer :: e2max=12
    integer :: lmax = -1
    logical :: test_mode = .false.
    logical :: count_memory = .false.
    logical :: ex

    ! input file name
    character(:), allocatable :: inputfile
    namelist /input/ zeta, basis, file_name, emax, e2max, lmax, test_mode, count_memory
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
    if(params%lmax==-1) params%lmax = params%emax

    params%count_memory = count_memory
  end subroutine InitInputParameters

  subroutine CopyInputParameters(params2, params1)
    class(InputParameters), intent(inout) :: params2
    type(InputParameters), intent(in) :: params1
    integer :: i


    params2%file_name= params1%file_name
    params2%emax     = params1%emax
    params2%lmax     = params1%lmax
    params2%e2max    = params1%e2max
    params2%test_mode= params1%test_mode
    params2%zeta     = params1%zeta
    params2%basis    = params1%basis
  end subroutine CopyInputParameters

  subroutine PrintInputParameters(params,unt)
    class(InputParameters), intent(in) :: params
    integer, intent(in), optional :: unt
    integer :: iunit, i
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
    write(iunit,'(a, i3, a, i3, a, i3)') '# emax: ', params%emax, ', e2max: ', params%e2max, &
        & ", lmax: ", params%lmax
    write(iunit,'(3a,f6.2)') '# basis: ', trim(params%basis), ", basis parameter: ", params%zeta
    write(iunit,*)
  end subroutine PrintInputParameters
end module AtHamilInput
