module EleSingleParticleState
  implicit none

  public :: EleSingleParticleOrbit
  public :: EleOrbits

  private :: SetEleSingleParticleOrbit
  private :: InitEleOrbits
  private :: FinEleOrbits
  private :: GetLabelFromIndex
  private :: GetIndexFromLabel

  type :: EleSingleParticleOrbit
    integer :: n = -1
    integer :: l = -1
    integer :: j = -1
    integer :: e = -1
    integer :: idx = -1
  contains
    procedure :: set => SetEleSingleParticleOrbit
  end type EleSingleParticleOrbit

  type :: EleOrbits
    integer, allocatable :: nlj2idx(:,:,:)
    type(EleSingleParticleOrbit), allocatable :: orb(:)
    logical :: is_constructed=.false.
    integer :: emax, lmax, norbs
    character(:), allocatable :: basis
  contains
    procedure :: init => InitEleOrbits
    procedure :: fin => FinEleOrbits
    procedure :: GetLabelFromIndex
    procedure :: GetIndexFromLabel
  end type EleOrbits

  character(1), private :: OrbitalAngMom(21)=['s',&
      & 'p','d','f','g','h','i','k','l','m','n',&
      & 'o','q','r','t','u','v','w','x','y','z']

contains

  subroutine FinEleOrbits(this)
    class(EleOrbits), intent(inout) :: this
    if(.not. this%is_constructed) return
#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In FinOrbitsIsospin'
#endif
    deallocate(this%orb)
    deallocate(this%nlj2idx)
    this%is_constructed = .false.
  end subroutine FinEleOrbits

  subroutine InitEleOrbits(this, emax, lmax, mode)
    class(EleOrbits), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    character(*), intent(in), optional :: mode
    integer :: i

#ifdef SingleParticleStateDebug
    write(*,'(a)') 'In InitEleOrbits:'
#endif

    this%basis = "HO"
    if(present(mode)) this%basis = mode

    select case(this%basis)
    case("HO", "ho", "LO-2nl", "lo-2nl")
      call set_orbits_ho(this, emax, lmax)
    case("AO", "ao")
      call set_orbits_ao(this, emax, lmax)
    case("LO", "lo")
      call set_orbits_lo(this, emax, lmax)
    case default
      write(*,"(2a)") "Unkown basis type: ", trim(this%basis)
      write(*,"(a)") "HO orbital is assumed"
      this%basis = "HO"
      call set_orbits_ho(this, emax, lmax)
    end select

    this%lmax = -1
    do i = 1, this%norbs
      this%lmax = max(this%lmax, this%orb(i)%l)
    end do
  end subroutine InitEleOrbits

  subroutine set_orbits_ho(this, emax, lmax)
    class(EleOrbits), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    integer :: e, l, n, s, j, cnt
    this%emax = emax
    this%lmax = emax
    if(present(lmax)) this%lmax = lmax
    allocate(this%nlj2idx(0:this%emax/2, 0:this%lmax, 1:2*this%lmax+1))
    this%nlj2idx(:,:,:) = 0
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
        end do
      end do
    end do
    this%norbs = cnt
    allocate(this%orb(this%norbs))
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        if(mod(e - l, 2) == 1) cycle
        n = (e - l) / 2
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
          call this%orb(cnt)%set(n,l,j,2*n+l,cnt)
          this%nlj2idx(n,l,j) = cnt
        end do
      end do
    end do
    this%is_constructed = .true.
  end subroutine set_orbits_ho

  subroutine set_orbits_ao(this, emax, lmax)
    class(EleOrbits), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    integer :: e, l, n, s, j, cnt
    this%emax = emax
    this%lmax = emax
    if(present(lmax)) this%lmax = lmax
    allocate(this%nlj2idx(0:this%emax, 0:this%lmax, 1:2*this%lmax+1))
    this%nlj2idx(:,:,:) = 0
    cnt = 0
    do e = 1, this%emax
      do l = 0, min(e-1,this%lmax)
        n = e
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
        end do
      end do
    end do
    this%norbs = cnt
    allocate(this%orb(this%norbs))
    cnt = 0
    do e = 1, this%emax
      do l = 0, min(e-1,this%lmax)
        n = e
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
          call this%orb(cnt)%set(n,l,j,e,cnt)
          this%nlj2idx(n,l,j) = cnt
        end do
      end do
    end do
    this%is_constructed = .true.
  end subroutine set_orbits_ao

  subroutine set_orbits_lo(this, emax, lmax)
    class(EleOrbits), intent(inout) :: this
    integer, intent(in) :: emax
    integer, intent(in), optional :: lmax
    integer :: e, l, n, s, j, cnt
    this%emax = emax
    this%lmax = emax
    if(present(lmax)) this%lmax = lmax
    allocate(this%nlj2idx(0:this%emax, 0:this%lmax, 1:2*this%lmax+1))
    this%nlj2idx(:,:,:) = 0
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        n = e - l
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
        end do
      end do
    end do
    this%norbs = cnt
    allocate(this%orb(this%norbs))
    cnt = 0
    do e = 0, this%emax
      do l = 0, min(e,this%lmax)
        n = e - l
        do s = -1, 1, 2
          j = 2*l + s
          if(j < 1) cycle
          cnt = cnt + 1
          call this%orb(cnt)%set(n,l,j,e,cnt)
          this%nlj2idx(n,l,j) = cnt
        end do
      end do
    end do
    this%is_constructed = .true.
  end subroutine set_orbits_lo

  subroutine SetEleSingleParticleOrbit(this, n, l, j, e, idx)
    class(EleSingleParticleOrbit), intent(inout) :: this
    integer, intent(in) :: n, l, j, e, idx
    this%n = n
    this%l = l
    this%j = j
    this%e = e
    this%idx = idx
#ifdef SingleParticleStateDebug
    write(*,'(a,i3,a,i3,a,i3,a,i3,a,i3)') &
        & 'index=', idx, ', n=', n, ', l=', l, ', j=', j, ', e=', e
#endif
  end subroutine SetEleSingleParticleOrbit

  function GetLabelFromIndex(this,idx) result(r)
    use AtLibrary, only: str
    class(EleOrbits), intent(in) :: this
    integer, intent(in) :: idx
    character(:), allocatable :: r
    integer :: n, l, j

    if(idx > this%norbs) then
      write(*,'(a)') "Warning in GetLabelFromIndexIsospin"
      r = ''
      return
    end if

    n = this%orb(idx)%n
    l = this%orb(idx)%l
    j = this%orb(idx)%j
    r = trim(str(n)) // OrbitalAngMom(l+1) // trim(str(j)) // '/2'
  end function GetLabelFromIndex

  function GetIndexFromLabel(this, label) result(r)
    class(EleOrbits), intent(in) :: this
    character(*), intent(in) :: label
    integer :: r
    character(:), allocatable :: str, str_n, str_l, str_j
    integer :: n, l, j, i
    r = 0

    if(label == '') then
      write(*,'(a)') "Warning in GetIndexFromLabelIsospin"
      return
    end if

    str = label(1:)
    str_n = ''
    do
      if(scan(str,'1234567890') /= 1) exit
      str_n = trim(str_n) // str(1:1)
      str = str(2:)
    end do
    str_l = str(1:1)

    str = str(2:)
    str_j = ''
    do
      if(scan(str,'1234567890') /= 1) exit
      str_j = trim(str_j) // str(1:1)
      str = str(2:)
    end do

    do i = 1, size(OrbitalAngMom)
      if(str_l == OrbitalAngMom(i)) then
        l = i - 1
        exit
      end if
      if(i == size(OrbitalAngMom)) then
        write(*,'(a)') "Error in GetIndexFromLabelIsospin: orbital angular momentum limitation"
        l = 0
        exit
      end if
    end do

    read(str_n,*) n
    read(str_j,*) j
    r = this%nlj2idx(n,l,j)
  end function GetIndexFromLabel
end module EleSingleParticleState

! main program for check
!program test
!  use SingleParticleState
!  type(Orbits) :: o
!  type(OrbitsIsospin) :: io
!  character(:), allocatable :: str
!  call o%init(4)
!  str = o%GetLabelFromIndex(1)
!  write(*,*) 1, o%GetIndexFromLabel(str)
!  str = o%GetLabelFromIndex(20)
!  write(*,*) 20, o%GetIndexFromLabel(str)
!  str = o%GetLabelFromIndex(15)
!  write(*,*) 15, o%GetIndexFromLabel(str)
!  call o%fin()
!
!
!  call o%init(4,2)
!  call o%fin()
!
!  call io%init(4)
!  str = io%GetLabelFromIndexIsospin(1)
!  write(*,*) 1, io%GetIndexFromLabelIsospin(str)
!  str = io%GetLabelFromIndexIsospin(20)
!  write(*,*) 20, io%GetIndexFromLabelIsospin(str)
!  str = io%GetLabelFromIndexIsospin(15)
!  write(*,*) 15, io%GetIndexFromLabelIsospin(str)
!  call io%fin()
!
!  call io%init(4,2)
!  call io%fin()
!end program test


