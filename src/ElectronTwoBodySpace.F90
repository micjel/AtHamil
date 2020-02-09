module ElectronTwoBodySpace
  use EleSingleParticleState
  implicit none
  type :: EleTwoBodyChan
    integer :: n_state
    integer :: j
    integer :: p
    integer, allocatable :: labels2n(:,:)
    integer, allocatable :: n2label1(:)
    integer, allocatable :: n2label2(:)
    integer, allocatable :: iphase(:,:)
  contains
    procedure :: init => InitEleTwoBodyChan
    procedure :: fin => FinEleTwoBodyChan
  end type EleTwoBodyChan

  type :: EleTwoBodySpace
    type(EleTwoBodyChan), allocatable :: jp(:)
    type(EleOrbits), pointer :: sps
    integer, allocatable :: jp2idx(:,:)
    integer :: NChan
    integer :: emax, e2max
    real(8) :: zeta ! length parameter is a = (a_0 / zeta)
    logical :: is_Constructed=.false.
  contains
    procedure :: init => InitEleTwoBodySpace
    procedure :: fin => FinEleTwoBodySpace
  end type EleTwoBodySpace

contains

  subroutine InitEleTwoBodySpace(this, sps, e2max, zeta)
    use AtLibrary, only: triag
    class(EleTwoBodySpace), intent(inout) :: this
    type(EleOrbits), target, intent(in) :: sps
    integer, intent(in) :: e2max
    real(8), intent(in) :: zeta
    integer :: ich
    integer, allocatable :: jj(:), pp(:), nn(:)
    this%emax = sps%emax
    this%e2max = e2max
    this%zeta = zeta
    this%sps => sps

    this%NChan = CountChannel(sps,e2max)
    allocate(jj(this%NChan))
    allocate(pp(this%NChan))
    allocate(nn(this%NChan))
    allocate(this%jp(this%NChan))
    allocate(this%jp2idx(0:2*sps%lmax+1,-1:1))
    this%jp2idx(:,:) = 0
    call SetChannelIndex(sps,e2max,jj,pp,nn)

    do ich = 1, this%NChan
      this%jp(ich)%n_state = nn(ich)
      this%jp2idx(jj(ich),pp(ich)) = ich
      call this%jp(ich)%Init(sps, jj(ich), pp(ich), this%e2max)
    end do
    deallocate(jj,pp,nn)
    this%is_Constructed = .true.
  end subroutine InitEleTwoBodySpace

  subroutine FinEleTwoBodySpace(this)
    class(EleTwoBodySpace), intent(inout) :: this
    integer :: ich
    if(.not. this%is_Constructed) return
    do ich = 1, this%NChan
      call this%jp(ich)%fin()
    end do
    deallocate(this%jp)
    deallocate(this%jp2idx)
    this%is_Constructed = .false.
    this%sps => null()
  end subroutine FinEleTwoBodySpace

  function CountChannel(sps, e2max) result(NChan)
    use AtLibrary, only: triag
    type(EleOrbits), target, intent(in) :: sps
    integer, intent(in) :: e2max
    type(EleSingleParticleOrbit), pointer :: o1, o2
    integer :: NChan
    integer :: ich, j, p, n
    integer :: i1, j1, l1
    integer :: i2, j2, l2
    ich = 0
    do j = 0, e2max + 1
      do p = 1, -1, -2
        n = 0
        do i1 = 1, sps%norbs
          o1 => sps%orb(i1)
          j1 = o1%j
          l1 = o1%l
          do i2 = 1, i1
            o2 => sps%orb(i2)
            j2 = o2%j
            l2 = o2%l
            if(o1%e + o2%e > e2max) cycle
            if(triag(j1, j2, 2*j)) cycle
            if((-1) ** (l1 + l2) /= p) cycle
            if(i1 == i2 .and. mod(j, 2) == 1) cycle
            n = n + 1
          end do
        end do
        if(n /= 0) ich = ich + 1
      end do
    end do
    NChan = ich
  end function CountChannel

  subroutine SetChannelIndex(sps, e2max,jj,pp,nn)
    use AtLibrary, only: triag
    type(EleOrbits), target, intent(in) :: sps
    type(EleSingleParticleOrbit), pointer :: o1, o2
    integer, intent(in) :: e2max
    integer, intent(inout) :: jj(:), pp(:), nn(:)
    integer :: ich, j, p, n
    integer :: i1, j1, l1
    integer :: i2, j2, l2
    ich = 0
    do j = 0, e2max + 1
      do p = 1, -1, -2
        n = 0
        do i1 = 1, sps%norbs
          o1 => sps%orb(i1)
          j1 = o1%j
          l1 = o1%l
          do i2 = 1, i1
            o2 => sps%orb(i2)
            j2 = o2%j
            l2 = o2%l
            if(o1%e + o2%e > e2max) cycle
            if(triag(j1, j2, 2*j)) cycle
            if((-1) ** (l1 + l2) /= p) cycle
            if(i1 == i2 .and. mod(j, 2) == 1) cycle
            n = n + 1
          end do
        end do
        if(n /= 0) then
          ich = ich + 1
          jj(ich) = j
          pp(ich) = p
          nn(ich) = n
        end if
      end do
    end do
  end subroutine SetChannelIndex

  subroutine InitEleTwoBodyChan(this, sps, j, p, e2max)
    use AtLibrary, only: triag
    class(EleTwoBodyChan), intent(inout) :: this
    type(EleOrbits), target, intent(in) :: sps
    type(EleSingleParticleOrbit), pointer :: o1, o2
    integer, intent(in) :: j, p, e2max
    integer :: n, m, i1, i2
    integer :: j1, l1, j2, l2

    n = this%n_state
    m = sps%norbs
    this%j = j
    this%p = p
    allocate(this%n2label1(n))
    allocate(this%n2label2(n))
    allocate(this%labels2n(m,m))
    allocate(this%iphase(m,m))
    this%n2label1 = 0
    this%n2label2 = 0
    this%labels2n = 0
    this%iphase = 0

    n = 0
    do i1 = 1, m
      o1 => sps%orb(i1)
      j1 = o1%j
      l1 = o1%l
      do i2 = 1, i1
        o2 => sps%orb(i2)
        j2 = o2%j
        l2 = o2%l
        if(o1%e + o2%e > e2max) cycle
        if(triag(j1, j2, 2*j)) cycle
        if((-1) ** (l1 + l2) /= p) cycle
        if(i1 == i2 .and. mod(j, 2) == 1) cycle
        n = n + 1
        this%n2label1(n) = i1
        this%n2label2(n) = i2
        this%labels2n(i1,i2) = n
        this%labels2n(i2,i1) = n
        this%iphase(i1,i2) = 1
        this%iphase(i2,i1) = -(-1) ** ((j1 + j2) / 2 - j)
      end do
    end do
  end subroutine InitEleTwoBodyChan

  subroutine FinEleTwoBodyChan(this)
    class(EleTwoBodyChan), intent(inout) :: this
    deallocate(this%n2label1)
    deallocate(this%n2label2)
    deallocate(this%labels2n)
    deallocate(this%iphase)
  end subroutine FinEleTwoBodyChan
end module ElectronTwoBodySpace
