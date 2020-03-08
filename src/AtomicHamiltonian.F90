module AtomicHamiltonian
  use omp_lib
  use LinAlgLib
  use ElectronTwoBodySpace
  use EleSingleParticleState
  implicit none

  type, extends(DMat) :: AtomicHamilChan
    logical :: Zero = .true.
  contains
  end type AtomicHamilChan

  type :: AtomicHamil
    type(AtomicHamilChan), allocatable :: ee_coulomb(:,:)
    type(AtomicHamilChan), allocatable :: ee_darwin(:,:)
    type(AtomicHamilChan), allocatable :: ee_spin_contact(:,:)
    type(AtomicHamilChan), allocatable :: ee_spin_orbit(:,:)
    type(AtomicHamilChan), allocatable :: ee_orbit_orbit(:,:)
    type(AtomicHamilChan), allocatable :: ee_spin_dipole(:,:)
    type(AtomicHamilChan) :: kinetic
    type(AtomicHamilChan) :: potential
    type(AtomicHamilChan) :: kinetic_p4  ! relativistic correction for kinetic term
    type(AtomicHamilChan) :: Darwin_term ! nucleus-electron contact interaction
    type(AtomicHamilChan) :: LS_term     ! nucleus-electron spin-orbit interaction
    type(AtomicHamilChan) :: S
    type(EleTwoBodySpace), pointer :: ms
    logical :: is_init = .false.
    integer :: jr = 0
    integer :: pr = 1
    character(256) :: OpName = "Hamil"
  contains
    procedure :: InitAtomicHamil
    procedure :: FinAtomicHamil
    procedure :: SetAtomicHamil
    procedure :: WriteAtomicHamil
    procedure :: GetFileNameAtomicHamil

    generic :: init => InitAtomicHamil
    generic :: fin => FinAtomicHamil
    generic :: set => SetAtomicHamil
    generic :: writef => WriteAtomicHamil
    generic :: GetFileName => GetFileNameAtomicHamil
  end type AtomicHamil

  real(8), private, allocatable :: rmesh(:), rwmesh(:)
  real(8), private, allocatable :: pmesh(:), pwmesh(:)
  integer, private :: NMesh = 100
  integer, private :: NMesh_Mom = 100
  real(8), private :: rmax = 0.d0
  real(8), private :: pmax = 0.d0
  real(8), private, allocatable :: rnl(:,:,:)
  real(8), private, allocatable :: rnl_mom(:,:,:)
  real(8), private :: a_0, Eh

  type, private :: FL
    real(8), allocatable :: F(:)
    real(8), allocatable :: Contact(:)
  end type FL

  type, private :: Fnl
    type(FL), allocatable :: bk(:,:)
    integer, allocatable :: idx(:,:,:,:)
    integer, allocatable :: n1(:)
    integer, allocatable :: n2(:)
    integer, allocatable :: l1(:)
    integer, allocatable :: l2(:)
    integer :: nidx
  contains
    procedure :: InitFNL
    procedure :: FinFNL
    procedure :: GetFNL
    procedure :: GetFNLContact
    generic :: init => InitFNL
    generic :: fin => FinFNL
  end type Fnl
  type(Fnl) :: Fintegral

contains

  subroutine FinFNL(this)
    class(Fnl), intent(inout) :: this
    integer :: bra, ket
    do bra = 1, this%nidx
      do ket = 1, this%nidx
        deallocate(this%bk(bra,ket)%F)
        deallocate(this%bk(bra,ket)%Contact)
      end do
    end do
    deallocate(this%bk)
    deallocate(this%idx)
    deallocate(this%n1)
    deallocate(this%n2)
    deallocate(this%l1)
    deallocate(this%l2)
  end subroutine FinFNL

  subroutine InitFNL(this, ms)
    use AtLibrary, only: gauss_legendre, laguerre_radial_wf_norm
    class(Fnl), intent(inout) :: this
    type(EleTwoBodySpace), intent(in) :: ms
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa
    integer :: nmin, nmax, lmin, lmax
    integer :: l1max, l2max
    integer :: a, L
    integer :: n1, l1, n2, l2
    integer :: n3, l3, n4, l4
    integer :: bra, ket
    integer :: cnt, i1, i2
    real(8) :: integral, r1, r2, log_r, zeta
    real(8), allocatable :: inter(:,:,:)
    real(8), allocatable, save :: r_1(:), w_1(:)
    real(8), allocatable, save :: r_2(:), w_2(:)
    real(8), allocatable, save :: Rnl1(:), Rnl2(:)
    real(8) :: rmax_
    integer :: n_1, n_2
    !$omp threadprivate(r_1, w_1, r_2, w_2, Rnl1, Rnl2)


    sps => ms%sps
    nmax = -1
    lmax = -1
    nmin = 1000
    lmin = 1000
    zeta = ms%zeta
    do a = 1, sps%norbs
      oa => sps%orb(a)
      nmin = min(nmin, oa%n)
      lmin = min(lmin, oa%l)
      nmax = max(nmax, oa%n)
      lmax = max(lmax, oa%l)
    end do

    allocate(this%idx(nmin:nmax,lmin:lmax,nmin:nmax,lmin:lmax))
    this%idx(:,:,:,:) = 0

    cnt = 0
    do n1 = nmin, nmax
      do n2 = nmin, nmax
        select case(sps%basis)
        case("AO","ao")
          l1max = min(n1 - 1, lmax)
          l2max = min(n2 - 1, lmax)
        case default
          l1max = lmax
          l2max = lmax
        end select
        do l1 = lmin, l1max
          do l2 = lmin, l2max
            cnt = cnt + 1
          end do
        end do
      end do
    end do
    this%nidx = cnt
    allocate(this%bk(this%nidx, this%nidx))
    allocate(this%n1(this%nidx))
    allocate(this%l1(this%nidx))
    allocate(this%n2(this%nidx))
    allocate(this%l2(this%nidx))
    cnt = 0
    do n1 = nmin, nmax
      do n2 = nmin, nmax
        select case(sps%basis)
        case("AO","ao")
          l1max = min(n1 - 1, lmax)
          l2max = min(n2 - 1, lmax)
        case default
          l1max = lmax
          l2max = lmax
        end select
        do l1 = lmin, l1max
          do l2 = lmin, l2max
            cnt = cnt + 1
            this%idx(n1,l1,n2,l2) = cnt
            this%n1(cnt) = n1
            this%l1(cnt) = l1
            this%n2(cnt) = n2
            this%l2(cnt) = l2
          end do
        end do
      end do
    end do

    write(*,"(a,i8)") "Storing e-e integrals n_{1}l_{1}n_{3}l_{3} combination: ", this%nidx

    do bra = 1, this%nidx
      n1 = this%n1(bra)
      l1 = this%l1(bra)
      n3 = this%n2(bra)
      l3 = this%l2(bra)
      do ket = 1, this%nidx
        n2 = this%n1(ket)
        l2 = this%l1(ket)
        n4 = this%n2(ket)
        l4 = this%l2(ket)

        lmin = max(abs(l1-l3)-1, abs(l2-l4)-1,0)
        lmax = min(abs(l1+l3), abs(l2+l4)) + 1
        allocate(this%bk(bra,ket)%F(lmin:lmax))
        allocate(this%bk(bra,ket)%Contact(lmin:lmax))
        this%bk(bra,ket)%F(:) = 0.d0
        this%bk(bra,ket)%Contact(:) = 0.d0
      end do
    end do
    !
    !return ! test
    !

    rmax_ = maxval(rmesh)
    allocate(inter(NMesh,0:2*lmax+1,this%nidx))
    inter(:,:,:) = 0.d0
    !$omp parallel
    !$omp do private(i2, r2, n_1, n_2, ket, n1, l1, n3, l3, i1, L, integral, r1, log_r)
    do i2 = 1, NMesh
      r2 = rmesh(i2) / zeta
      n_1 = max(int(rmesh(i2) / rmax_ * dble(NMesh)),NMesh/20)
      n_2 = NMesh - n_1
      call gauss_legendre(0.d0, rmesh(i2), r_1, w_1, n_1)
      call gauss_legendre(rmesh(i2), rmax_, r_2, w_2, n_2)
      allocate( Rnl1(n_1) )
      allocate( Rnl2(n_2) )
      do ket = 1, this%nidx
        n1 = this%n1(ket)
        l1 = this%l1(ket)
        n3 = this%n2(ket)
        l3 = this%l2(ket)
        do i1 = 1, n_1
          Rnl1(i1) = laguerre_radial_wf_norm(n1,l1,1.d0,r_1(i1)) * laguerre_radial_wf_norm(n3,l3,1.d0,r_1(i1)) * w_1(i1)
        end do
        do i1 = 1, n_2
          Rnl2(i1) = laguerre_radial_wf_norm(n1,l1,1.d0,r_2(i1)) * laguerre_radial_wf_norm(n3,l3,1.d0,r_2(i1)) * w_2(i1)
        end do

        do L = abs(l1-l3)-1, l1+l3+1
          if(L < 0) cycle
          if(mod(l1 + l3 + L, 2) == 1) cycle

          integral = 0.d0
          do i1 = 1, n_1
            r1 = r_1(i1) / zeta
            log_r = dble(L) * log(r1) - dble(L+1) * log(r2)
            integral = integral + exp(log_r) * Rnl1(i1)
          end do
          do i1 = 1, n_2
            r1 = r_2(i1) / zeta
            log_r = dble(L) * log(r2) - dble(L+1) * log(r1)
            integral = integral + exp(log_r) * Rnl2(i1)
          end do
          inter(i2,L,ket) = integral
        end do
      end do
      deallocate(r_1, w_1, r_2, w_2, Rnl1, Rnl2)
    end do
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do private(bra, n1, l1, n3, l3, ket, n2, l2, n4, l4, lmin, lmax, &
    !$omp &          L, integral, i2)
    do bra = 1, this%nidx
      n1 = this%n1(bra)
      l1 = this%l1(bra)
      n3 = this%n2(bra)
      l3 = this%l2(bra)
      do ket = 1, this%nidx
        n2 = this%n1(ket)
        l2 = this%l1(ket)
        n4 = this%n2(ket)
        l4 = this%l2(ket)
        lmin = max(abs(l1-l3)-1, abs(l2-l4)-1, 0)
        lmax = min(l1+l3, l2+l4) + 1
        do L = lmin, lmax
          if(L < 0) cycle
          if(mod(l1 + l3 + L, 2) == 1) cycle
          if(mod(l2 + l4 + L, 2) == 1) cycle
          integral = 0.d0
          do i2 = 1, NMesh
            integral = integral + inter(i2,L,ket) * rwmesh(i2) * &
                & rnl(i2,n1,l1) * rnl(i2,n3,l3)
          end do
          this%bk(bra,ket)%F(L) = integral
          !write(*,"(9i3,f12.6)") n1,l1,n2,l2,n3,l3,n4,l4,L,integral
          !write(*,*) bra, ket, L, integral
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    !$omp parallel
    !$omp do private(bra, n1, l1, n3, l3, ket, n2, l2, n4, l4, lmin, lmax, &
    !$omp &          L, integral, i1)
    do bra = 1, this%nidx
      n1 = this%n1(bra)
      l1 = this%l1(bra)
      n3 = this%n2(bra)
      l3 = this%l2(bra)
      do ket = 1, this%nidx
        n2 = this%n1(ket)
        l2 = this%l1(ket)
        n4 = this%n2(ket)
        l4 = this%l2(ket)
        lmin = max(abs(l1-l3)-1, abs(l2-l4)-1, 0)
        lmax = min(l1+l3, l2+l4) + 1
        do L = lmin, lmax
          if(L < 0) cycle
          if(mod(l1 + l3 + L, 2) == 1) cycle
          if(mod(l2 + l4 + L, 2) == 1) cycle
          integral = 0.d0
          do i1 = 1, NMesh
            integral = integral + rwmesh(i1) * &
                & rnl(i1,n1,l1) * rnl(i1,n3,l3) * rnl(i1,n2,l2) * rnl(i1,n4,l4) * (zeta/rmesh(i1))**2
          end do
          this%bk(bra,ket)%Contact(L) = integral
          !!$omp critical
          !write(*,"(9i3,f12.6)") n1,l1,n2,l2,n3,l3,n4,l4,L,integral
          !!$omp end critical
          !write(*,*) bra, ket, L, integral
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(inter)
  end subroutine InitFNL

  function GetFNL(this, n1, l1, n2, l2, n3, l3, n4, l4, L) result(r)
    class(Fnl), intent(in) :: this
    integer, intent(in) :: n1, l1, n2, l2, n3, l3, n4, l4, L
    integer :: bra, ket
    real(8) :: r
    bra = this%idx(n1,l1,n3,l3)
    ket = this%idx(n2,l2,n4,l4)
    r = this%bk(bra,ket)%F(L)
  end function GetFNL

  function GetFNLContact(this, n1, l1, n2, l2, n3, l3, n4, l4, L) result(r)
    class(Fnl), intent(in) :: this
    integer, intent(in) :: n1, l1, n2, l2, n3, l3, n4, l4, L
    integer :: bra, ket
    real(8) :: r
    bra = this%idx(n1,l1,n3,l3)
    ket = this%idx(n2,l2,n4,l4)
    r = this%bk(bra,ket)%Contact(L)
  end function GetFNLContact

  subroutine FinAtomicHamil(this)
    class(AtomicHamil), intent(inout) :: this
    integer :: chbra, chket
    call Fintegral%FinFNL()
    do chbra = 1, this%ms%NChan
      do chket = 1, this%ms%NChan
        call this%ee_coulomb(chbra,chket)%fin()
        call this%ee_darwin(chbra,chket)%fin()
        call this%ee_spin_contact(chbra,chket)%fin()
        call this%ee_spin_orbit(chbra,chket)%fin()
        call this%ee_orbit_orbit(chbra,chket)%fin()
        call this%ee_spin_dipole(chbra,chket)%fin()
      end do
    end do
    deallocate(this%ee_coulomb)
    deallocate(this%ee_darwin)
    deallocate(this%ee_spin_contact)
    deallocate(this%ee_spin_orbit)
    deallocate(this%ee_orbit_orbit)
    deallocate(this%ee_spin_dipole)
    call this%kinetic%fin()
    call this%potential%fin()
    call this%kinetic_p4%fin()
    call this%Darwin_term%fin()
    call this%LS_term%fin()
    this%ms => null()
    this%is_init = .false.
    deallocate(rmesh)
    deallocate(rwmesh)
    deallocate(pmesh)
    deallocate(pwmesh)
    deallocate(rnl)
    deallocate(rnl_mom)
  end subroutine FinAtomicHamil

  subroutine InitAtomicHamil(this, two)
    use AtLibrary, only: m_e, hc, alpha
    class(AtomicHamil), intent(inout) :: this
    type(EleTwoBodySpace), target, intent(in) :: two
    integer :: chbra, chket
    integer :: nb, jb, pb
    integer :: nk, jk, pk
    if(allocated(this%ee_coulomb)) call this%fin()
    this%ms => two
    allocate(this%ee_coulomb(two%NChan, two%NChan))
    allocate(this%ee_darwin(two%NChan, two%NChan))
    allocate(this%ee_spin_contact(two%NChan, two%NChan))
    allocate(this%ee_spin_orbit(two%NChan, two%NChan))
    allocate(this%ee_orbit_orbit(two%NChan, two%NChan))
    allocate(this%ee_spin_dipole(two%NChan, two%NChan))
    do chbra = 1, two%NChan
      nb = this%ms%jp(chbra)%n_state
      jb = this%ms%jp(chbra)%j
      pb = this%ms%jp(chbra)%p
      do chket = 1, two%NChan
        nk = this%ms%jp(chket)%n_state
        jk = this%ms%jp(chket)%j
        pk = this%ms%jp(chket)%p
        if(jb /= jk) cycle
        if(pb /= pk) cycle
        this%ee_coulomb(chbra,chket)%Zero = .false.
        this%ee_darwin(chbra,chket)%Zero = .false.
        this%ee_spin_contact(chbra,chket)%Zero = .false.
        this%ee_spin_orbit(chbra,chket)%Zero = .false.
        this%ee_orbit_orbit(chbra,chket)%Zero = .false.
        this%ee_spin_dipole(chbra,chket)%Zero = .false.
        call this%ee_coulomb(chbra, chket)%zeros(nb,nk)
        call this%ee_darwin(chbra, chket)%zeros(nb,nk)
        call this%ee_spin_contact(chbra, chket)%zeros(nb,nk)
        call this%ee_spin_orbit(chbra, chket)%zeros(nb,nk)
        call this%ee_orbit_orbit(chbra, chket)%zeros(nb,nk)
        call this%ee_spin_dipole(chbra, chket)%zeros(nb,nk)
      end do
    end do
    this%kinetic%Zero = .false.
    this%potential%Zero = .false.
    call this%kinetic%zeros(two%sps%norbs, two%sps%norbs)
    call this%potential%zeros(two%sps%norbs, two%sps%norbs)
    call this%kinetic_p4%zeros(two%sps%norbs, two%sps%norbs)
    call this%Darwin_term%zeros(two%sps%norbs, two%sps%norbs)
    call this%LS_term%zeros(two%sps%norbs, two%sps%norbs)
    call this%S%eye(two%sps%norbs)
    this%is_init = .true.
    a_0 = hc / (m_e * 1.d3) * alpha ! bohr radius unit of nm nuclear mass infinity limit
    Eh = hc / (a_0 * alpha)         ! hartree
    select case(two%sps%basis)
    case("STO", "sto")
    case("GTO", "gto")
    end select

  end subroutine InitAtomicHamil

  function GetFileNameAtomicHamil(this) result(f)
    use AtLibrary, only: str
    class(AtomicHamil), intent(inout) :: this
    type(EleTwoBodySpace), pointer :: ms
    type(EleOrbits), pointer :: sps
    character(:), allocatable :: f
    ms => this%ms
    sps => ms%sps
    f = "AtomicHamil_" // trim(sps%basis) // "_zeta" // trim(str(ms%zeta))
    f = trim(f) // "_emax" // trim(str(sps%emax))
    f = trim(f) // "_e2max" // trim(str(ms%e2max)) // "_lmax" // trim(str(sps%lmax))
    f = trim(f) // ".snt"
  end function GetFileNameAtomicHamil

  subroutine norm_check(wf, wmesh, basis)
    real(8), intent(in) :: wf(:,:,:)
    real(8), intent(in) :: wmesh(:)
    character(*), intent(in) :: basis
    real(8) :: ovlp
    integer :: nbra, nket, l, i, lmax

    lmax = size(wf,3)
    do nbra = 1, size(wf,2)
      if(basis == "AO") lmax = min(nbra-1,size(wf,3))
      do nket = 1, size(wf,2)
        if(basis == "AO") lmax = min(nket-1,size(wf,3))
        if(nbra > nket) cycle
        do l = 1, lmax
          if(basis == "HO" .and. 2*nbra+l > 2*(size(wf,2)-1)) cycle
          if(basis == "HO" .and. 2*nket+l > 2*(size(wf,2)-1)) cycle

          ovlp = 0.d0
          do i = 1, size(wmesh)
            ovlp = ovlp + wmesh(i) * wf(i,nbra,l) * wf(i,nket,l)
          end do
          if(nbra == nket .and. abs(1.d0-ovlp) > 1.d-4) then
            write(*,"(a,3i3,f12.8)") "Warning: radial wave function norm", nbra,nket,l,ovlp
          end if

          if(nbra /= nket .and. abs(ovlp) > 1.d-4) then
            write(*,"(a,3i3,f12.8)") "Warning: radial wave function norm", nbra,nket,l,ovlp
          end if
        end do
      end do
    end do
  end subroutine norm_check

  subroutine SetAtomicHamil(this, NMesh_in, rmax_in, NMesh_mom_in, pmax_in)
    use AtLibrary, only: gauss_legendre, &
        & ho_radial_wf_norm, hydrogen_radial_wf_norm, &
        & hydrogen_radial_wf_mom_norm, &
        & laguerre_radial_wf_norm, &
        & fixed_point_quadrature, ln_gamma, laguerre, &
        & Mom_laguerre_radial_wf_norm
    class(AtomicHamil), intent(inout) :: this
    integer, intent(in) :: NMesh_in, NMesh_Mom_in
    real(8), intent(in) :: rmax_in, pmax_in
    type(EleTwoBodySpace), pointer :: two
    real(8) :: ti
    integer :: n, l, i

    NMesh = NMesh_in
    NMesh_mom = NMesh_mom_in
    rmax = rmax_in
    pmax = pmax_in

    two => this%ms
    ti = omp_get_wtime()
    select case(two%sps%basis)
    case("HO", "ho")
      write(*,*) "Using Gauss-Legendre quadrature"
      call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
      allocate(rnl(NMesh, 0:two%sps%emax/2, 0:two%sps%emax))
      rnl(:,:,:) = 0.d0
      do n = 0, two%sps%emax/2
        do l = 0, two%sps%emax
          do i = 1, NMesh
            rnl(i,n,l) = ho_radial_wf_norm(n, l, 1.d0, rmesh(i))
          end do
        end do
      end do

      call norm_check(rnl, rwmesh, "HO")
    case("AO", "ao")
      write(*,*) "Using Gauss-Legendre quadrature"
      call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
      call gauss_legendre(0.d0, pmax, pmesh, pwmesh, NMesh_Mom)
      allocate(rnl(NMesh, 1:two%sps%emax, 0:two%sps%lmax))
      allocate(rnl_mom(NMesh_Mom, 1:two%sps%emax, 0:two%sps%lmax))
      rnl(:,:,:) = 0.d0
      rnl_mom(:,:,:) = 0.d0
      do n = 1, two%sps%emax
        do l = 0, min(n-1,two%sps%lmax)
          do i = 1, NMesh
            rnl(i,n,l) = hydrogen_radial_wf_norm(n,l,1.d0,rmesh(i))
          end do
          do i = 1, NMesh_Mom
            rnl_mom(i,n,l) = hydrogen_radial_wf_mom_norm(n,l,1.d0,pmesh(i))
          end do
        end do
      end do
      call norm_check(rnl, rwmesh, "AO")
      call norm_check(rnl_mom, pwmesh, "AO")
    case("LO", "lo")
#ifdef gauss_laguerre
      write(*,*) "Using Gauss-Laguerre quadrature"
      call fixed_point_quadrature("laguerre", NMesh, rmesh, rwmesh, weight_renorm=.false., &
        & a_in=0.d0, b_in=2.d0, alpha_in=0.d0)
#else
      write(*,*) "Using Gauss-Legendre quadrature"
      call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
#endif
      call gauss_legendre(0.d0, pmax, pmesh, pwmesh, NMesh_Mom)
      allocate(rnl(NMesh, 0:two%sps%emax, 0:two%sps%lmax))
      allocate(rnl_mom(NMesh_Mom, 0:two%sps%emax, 0:two%sps%lmax))
      rnl(:,:,:) = 0.d0
      do n = 0, two%sps%emax
        do l = 0, two%sps%lmax
          do i = 1, NMesh
#ifdef gauss_laguerre
            rnl(i,n,l) = exp( 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+2*l+3))) * &
              & laguerre(n,dble(2*l+2),2.d0*rmesh(i)) * (2.d0*rmesh(i))**(l+1) * sqrt(2.d0)
#else
            rnl(i,n,l) = laguerre_radial_wf_norm(n, l, 1.d0, rmesh(i))
#endif
          end do
          do i = 1, NMesh_Mom
            rnl_mom(i,n,l) = Mom_laguerre_radial_wf_norm(n, l, 1.d0, pmesh(i))
          end do
        end do
      end do
      call norm_check(rnl, rwmesh, "lo")
      call norm_check(rnl_mom, pwmesh, "lo")
    case("LO-2nl", "lo-2nl")
#ifdef gauss_laguerre
      write(*,*) "Using Gauss-Laguerre quadrature"
      call fixed_point_quadrature("laguerre", NMesh, rmesh, rwmesh, weight_renorm=.false., &
          & a_in=0.d0, b_in=2.d0, alpha_in=0.d0)
#else
      write(*,*) "Using Gauss-Legendre quadrature"
      call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
#endif
      call gauss_legendre(0.d0, pmax, pmesh, pwmesh, NMesh_Mom)
      allocate(rnl(NMesh, 0:two%sps%emax, 0:two%sps%lmax))
      allocate(rnl_mom(NMesh_Mom, 1:two%sps%emax, 0:two%sps%lmax))
      rnl(:,:,:) = 0.d0
      do n = 0, two%sps%emax/2
        do l = 0, two%sps%lmax
          do i = 1, NMesh
#ifdef gauss_laguerre
            rnl(i,n,l) = exp( 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+2*l+3))) * &
              & laguerre(n,dble(2*l+2),2.d0*rmesh(i)) * (2.d0*rmesh(i))**(l+1) * sqrt(2.d0)
#else
            rnl(i,n,l) = laguerre_radial_wf_norm(n, l, 1.d0, rmesh(i))
#endif
          end do
          do i = 1, NMesh_Mom
            rnl_mom(i,n,l) = Mom_laguerre_radial_wf_norm(n, l, 1.d0, pmesh(i))
          end do
        end do
      end do
      call norm_check(rnl, rwmesh, "HO")
      call norm_check(rnl_mom, pwmesh, "HO")
    end select
    write(*,"(a,f12.6,a)") "Basis storing: ", omp_get_wtime() - ti, " sec"
    ti = omp_get_wtime()
    select case(two%sps%basis)
    case("HO","ho","AO","ao","LO","lo","LO-2nl")
      call Fintegral%init(two)
    end select
    write(*,"(a,f12.6,a)") "Radial integral storing:  ", omp_get_wtime() - ti, " sec"

    ti = omp_get_wtime()
    call set_one_body_part(this)
    call set_two_body_part(this)
    write(*,"(a,f12.6,a)") "SetHamil:  ", omp_get_wtime() - ti, " sec"
  end subroutine SetAtomicHamil

  subroutine set_one_body_part(this)
    use AtLibrary, only: gauss_legendre, ho_radial_wf_norm, hc, m_e
    class(AtomicHamil), intent(inout) :: this
    integer :: a, b, i
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob
    type(DMat) :: kine, U
    real(8) :: r, hw, zeta

    zeta = this%ms%zeta
    hw = hc **2 * zeta**2 / m_e*1.d3 ! in a.u.
    sps => this%ms%sps
    call kine%zeros(sps%norbs, sps%norbs)
    call U%zeros(sps%norbs, sps%norbs)
    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        r = 0.d0
        select case(sps%basis)
        case("HO", "ho")
          if(abs(oa%n - ob%n) > 1) cycle
          if(oa%n == ob%n) r = dble(2 * oa%n + oa%l) + 1.5d0
          if(oa%n == ob%n-1) r = dsqrt(dble(oa%n + 1) * (dble(oa%n + oa%l) + 1.5d0))
          if(oa%n == ob%n+1) r = dsqrt(dble(ob%n + 1) * (dble(ob%n + ob%l) + 1.5d0))
          r = r * hw * 0.5d0
        case("AO","ao")
          do i = 1, NMesh_Mom
            r = r + rnl_mom(i,oa%n,oa%l) * rnl_mom(i,ob%n,ob%l) * &
                & pwmesh(i) * (pmesh(i)*zeta)**2 * 0.5d0
          end do
        case("LO","lo","LO-2nl","lo-2nl")
          r = kinetic_energy_laguerre_func(oa%n,ob%n,oa%l)
          !r = r * hw / Eh
          r = r * this%ms%zeta**2
        end select
        kine%m(a,b) = r
        kine%m(b,a) = r
      end do
    end do

    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        r = 0.d0
        do i = 1, NMesh
          r = r + rnl(i,oa%n,oa%l) * rnl(i,ob%n,ob%l) * rwmesh(i) * zeta/ rmesh(i)
        end do
        U%m(a, b) = r
        U%m(b, a) = r
      end do
    end do

    this%kinetic%DMat = kine
    this%potential%DMat = U
    call set_kinetic_relativistic_correction(this)
    call set_Darwin_term(this)
    call set_LS_term(this)
    call kine%fin()
    call U%fin()
  end subroutine set_one_body_part

  subroutine set_kinetic_relativistic_correction(this)
    use AtLibrary, only: alpha
    class(AtomicHamil), intent(inout) :: this
    integer :: a, b, i
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob
    real(8) :: r, zeta

    zeta = this%ms%zeta
    sps => this%ms%sps
    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        r = 0.d0
        select case(sps%basis)
        case("LO","lo","LO-2nl","lo-2nl")
          do i = 1, NMesh_Mom
            r = r + pwmesh(i) * rnl_mom(i,oa%n,oa%l) * rnl_mom(i,ob%n,ob%l) * &
              & (pmesh(i)*zeta)**4 * 0.125d0 / alpha**2
          end do
        case default
          write(*,*) "Relativistic correction with the selected basis function is not implemented"
          return
        end select
        this%kinetic_p4%m(a,b) = -r
        this%kinetic_p4%m(b,a) = -r
      end do
    end do
  end subroutine set_kinetic_relativistic_correction

  subroutine set_Darwin_term(this)
    !
    ! Assuming the nucleus is a point charge
    !
    use AtLibrary, only: pi, laguerre_radial_wf, alpha
    class(AtomicHamil), intent(inout) :: this
    integer :: a, b
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob
    real(8) :: r, zeta

    zeta = this%ms%zeta
    sps => this%ms%sps
    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        if(oa%l /= 0) cycle
        r = 0.d0
        select case(sps%basis)
        case("LO","lo","LO-2nl","lo-2nl")
          r = 0.125d0 * laguerre_radial_wf(oa%n,0,1.d0/zeta,0.d0) * laguerre_radial_wf(ob%n,0,1.d0/zeta,0.d0) / alpha**2
        case default
          write(*,*) "Darwin term with the selected basis function is not implemented"
          return
        end select
        this%Darwin_term%m(a,b) = r
        this%Darwin_term%m(b,a) = r
      end do
    end do
  end subroutine set_Darwin_term

  subroutine set_LS_term(this)
    use AtLibrary, only: g_s, alpha
    class(AtomicHamil), intent(inout) :: this
    integer :: a, b, i
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob
    real(8) :: r, zeta, ls

    zeta = this%ms%zeta
    sps => this%ms%sps
    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        if(oa%l == 0) cycle
        ls = dble(oa%j)*0.5d0*(dble(oa%j)*0.5d0+1.d0) - dble(oa%l*(oa%l+1)) - 0.75d0
        r = 0.d0
        select case(sps%basis)
        case("LO","lo","LO-2nl","lo-2nl")
          do i = 1, NMesh
            r = r + rnl(i,oa%n,oa%l) * rnl(i,ob%n,ob%l) * rwmesh(i) * zeta**3/ rmesh(i)**3
          end do
          r = r * ls * 0.25d0 * g_s / alpha**2
        case default
          write(*,*) "Darwin term with the selected basis function is not implemented"
          return
        end select
        this%LS_term%m(a,b) = r
        this%LS_term%m(b,a) = r
      end do
    end do
  end subroutine set_LS_term

  subroutine set_two_body_part(this)
    class(AtomicHamil), intent(inout) :: this

    call set_ee_coulomb(this)
    call set_ee_darwin(this)
  end subroutine set_two_body_part

  subroutine set_ee_coulomb(this)
    class(AtomicHamil), intent(inout) :: this
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: bra, ket, a, b, c, d, J
    integer :: ch
    real(8) :: rabcd, rbacd, rabdc, rbadc, norm

    sps => this%ms%sps
    do ch = 1, this%ms%NChan
      J = this%ms%jp(ch)%j
      !$omp parallel
      !$omp do private(bra, a, b, oa, ob, ket, c, d, oc, od, norm, &
      !$omp &          rabcd, rbacd, rabdc, rbadc)
      do bra = 1, this%ms%jp(ch)%n_state
        a = this%ms%jp(ch)%n2label1(bra)
        b = this%ms%jp(ch)%n2label2(bra)
        oa => sps%orb(a)
        ob => sps%orb(b)
        do ket = 1, bra
          c = this%ms%jp(ch)%n2label1(ket)
          d = this%ms%jp(ch)%n2label2(ket)
          oc => sps%orb(c)
          od => sps%orb(d)

          norm = 1.d0
          if(a == b) norm = norm * dsqrt(0.5d0)
          if(c == d) norm = norm * dsqrt(0.5d0)
          rabcd = ee_interaction(oa, ob, oc, od, J)
          rbacd = ee_interaction(ob, oa, oc, od, J) * (-1.d0)**((oa%j+ob%j)/2-J-1)
          rabdc = ee_interaction(oa, ob, od, oc, J) * (-1.d0)**((oc%j+od%j)/2-J-1)
          rbadc = ee_interaction(ob, oa, od, oc, J) * (-1.d0)**((oa%j+ob%j+oc%j+od%j)/2)
          this%ee_coulomb(ch,ch)%m(bra,ket) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
          this%ee_coulomb(ch,ch)%m(ket,bra) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
  end subroutine set_ee_coulomb

  subroutine set_ee_darwin(this)
    use AtLibrary, only: alpha
    class(AtomicHamil), intent(inout) :: this
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: bra, ket, a, b, c, d, J
    integer :: ch
    real(8) :: rabcd, rbacd, rabdc, rbadc, norm

    sps => this%ms%sps
    do ch = 1, this%ms%NChan
      J = this%ms%jp(ch)%j
      !$omp parallel
      !$omp do private(bra, a, b, oa, ob, ket, c, d, oc, od, norm, &
      !$omp &          rabcd, rbacd, rabdc, rbadc)
      do bra = 1, this%ms%jp(ch)%n_state
        a = this%ms%jp(ch)%n2label1(bra)
        b = this%ms%jp(ch)%n2label2(bra)
        oa => sps%orb(a)
        ob => sps%orb(b)
        do ket = 1, bra
          c = this%ms%jp(ch)%n2label1(ket)
          d = this%ms%jp(ch)%n2label2(ket)
          oc => sps%orb(c)
          od => sps%orb(d)

          norm = 1.d0
          if(a == b) norm = norm * dsqrt(0.5d0)
          if(c == d) norm = norm * dsqrt(0.5d0)
          rabcd = ee_interaction_darwin(oa, ob, oc, od, J)
          rbacd = ee_interaction_darwin(ob, oa, oc, od, J) * (-1.d0)**((oa%j+ob%j)/2-J-1)
          rabdc = ee_interaction_darwin(oa, ob, od, oc, J) * (-1.d0)**((oc%j+od%j)/2-J-1)
          rbadc = ee_interaction_darwin(ob, oa, od, oc, J) * (-1.d0)**((oa%j+ob%j+oc%j+od%j)/2)
          this%ee_darwin(ch,ch)%m(bra,ket) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
          this%ee_darwin(ch,ch)%m(ket,bra) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
        end do
      end do
      !$omp end do
      !$omp end parallel
      this%ee_darwin(ch,ch)%m(:,:) = this%ee_darwin(ch,ch)%m(:,:) * 0.25d0 / alpha**2
    end do
  end subroutine set_ee_darwin

  subroutine set_ee_spin_contact(this)
    use AtLibrary, only: alpha
    class(AtomicHamil), intent(inout) :: this
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: bra, ket, a, b, c, d, J
    integer :: ch
    real(8) :: rabcd, rbacd, rabdc, rbadc, norm

    sps => this%ms%sps
    do ch = 1, this%ms%NChan
      J = this%ms%jp(ch)%j
      !$omp parallel
      !$omp do private(bra, a, b, oa, ob, ket, c, d, oc, od, norm, &
      !$omp &          rabcd, rbacd, rabdc, rbadc)
      do bra = 1, this%ms%jp(ch)%n_state
        a = this%ms%jp(ch)%n2label1(bra)
        b = this%ms%jp(ch)%n2label2(bra)
        oa => sps%orb(a)
        ob => sps%orb(b)
        do ket = 1, bra
          c = this%ms%jp(ch)%n2label1(ket)
          d = this%ms%jp(ch)%n2label2(ket)
          oc => sps%orb(c)
          od => sps%orb(d)

          norm = 1.d0
          if(a == b) norm = norm * dsqrt(0.5d0)
          if(c == d) norm = norm * dsqrt(0.5d0)
          rabcd = ee_interaction_spin_contact(oa, ob, oc, od, J)
          rbacd = ee_interaction_spin_contact(ob, oa, oc, od, J) * (-1.d0)**((oa%j+ob%j)/2-J-1)
          rabdc = ee_interaction_spin_contact(oa, ob, od, oc, J) * (-1.d0)**((oc%j+od%j)/2-J-1)
          rbadc = ee_interaction_spin_contact(ob, oa, od, oc, J) * (-1.d0)**((oa%j+ob%j+oc%j+od%j)/2)
          this%ee_spin_contact(ch,ch)%m(bra,ket) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
          this%ee_spin_contact(ch,ch)%m(ket,bra) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
        end do
      end do
      !$omp end do
      !$omp end parallel
      this%ee_spin_contact(ch,ch)%m(:,:) = this%ee_spin_contact(ch,ch)%m(:,:) * 2.d0 / (3.d0*alpha**2)
    end do
  end subroutine set_ee_spin_contact

  function ee_interaction(oa, ob, oc, od, J) result(r)
    use AtLibrary, only: tjs, sjs
    type(EleSingleParticleOrbit), intent(in) :: oa, ob, oc, od
    integer, intent(in) :: J
    real(8) :: r
    integer :: Lmin, Lmax, L
    real(8) :: integral
    Lmin = max(abs(oa%j-oc%j), abs(ob%j-od%j))/2
    Lmax = min(   (oa%j+oc%j),    (ob%j+od%j))/2
    r = 0.d0
    do L = Lmin, Lmax
      if(mod(oa%l + oc%l + L, 2) == 1) cycle
      if(mod(ob%l + od%l + L, 2) == 1) cycle

      integral = Fintegral%GetFNL(oa%n,oa%l,ob%n,ob%l,oc%n,oc%l,od%n,od%l,L)

      r = r + integral * sjs(oa%j, ob%j, 2*J, od%j, oc%j, 2*L) * &
          & tjs(oa%j, 2*L, oc%j, -1, 0, 1) * &
          & tjs(ob%j, 2*L, od%j, -1, 0, 1)
    end do
    r = r * dsqrt(dble(oa%j+1) *dble(ob%j+1) * dble(oc%j+1) * dble(od%j+1)) * &
        &   (-1.d0) ** ((oa%j + oc%j) /2 + J)
  end function ee_interaction

  function ee_interaction_darwin(oa, ob, oc, od, J) result(r)
    use AtLibrary, only: tjs, sjs
    type(EleSingleParticleOrbit), intent(in) :: oa, ob, oc, od
    integer, intent(in) :: J
    real(8) :: r
    integer :: Lmin, Lmax, L
    real(8) :: integral

    Lmin = max(abs(oa%j-oc%j), abs(ob%j-od%j))/2
    Lmax = min(   (oa%j+oc%j),    (ob%j+od%j))/2
    r = 0.d0
    do L = Lmin, Lmax
      if(mod(oa%l + oc%l + L, 2) == 1) cycle
      if(mod(ob%l + od%l + L, 2) == 1) cycle
      integral = Fintegral%GetFNLContact(oa%n,oa%l,ob%n,ob%l,oc%n,oc%l,od%n,od%l,L)

      r = r + integral * sjs(oa%j, ob%j, 2*J, od%j, oc%j, 2*L) * &
          & tjs(oa%j, 2*L, oc%j, -1, 0, 1) * &
          & tjs(ob%j, 2*L, od%j, -1, 0, 1)
    end do
    r = r * dsqrt(dble(oa%j+1) *dble(ob%j+1) * dble(oc%j+1) * dble(od%j+1)) * &
        &   (-1.d0) ** ((oa%j + oc%j) /2 + J)
  end function ee_interaction_darwin

  function ee_interaction_spin_contact(oa, ob, oc, od, J) result(r)
    use AtLibrary, only: tjs, sjs
    type(EleSingleParticleOrbit), intent(in) :: oa, ob, oc, od
    integer, intent(in) :: J
    real(8) :: r, fk
    integer :: Lmin, Lmax, L, k
    real(8) :: integral

    Lmin = max(abs(oa%j-oc%j), abs(ob%j-od%j))/2
    Lmax = min(   (oa%j+oc%j),    (ob%j+od%j))/2
    r = 0.d0
    do L = Lmin, Lmax
      if(mod(oa%l + oc%l + L, 2) == 1) cycle
      if(mod(ob%l + od%l + L, 2) == 1) cycle
      !fk = 0.d0
      !do k = abs(L-1), L+1
      !
      !end do
      !integral = Fintegral%GetFNLContact(oa%n,oa%l,ob%n,ob%l,oc%n,oc%l,od%n,od%l,L)

      !r = r + integral * sjs(oa%j, ob%j, 2*J, od%j, oc%j, 2*L) * &
      !    & tjs(oa%j, 2*L, oc%j, -1, 0, 1) * &
      !    & tjs(ob%j, 2*L, od%j, -1, 0, 1)
    end do
    r = r * dsqrt(dble(oa%j+1) *dble(ob%j+1) * dble(oc%j+1) * dble(od%j+1)) * &
        &   (-1.d0) ** ((oa%j + oc%j) /2 + J)
  end function ee_interaction_spin_contact

  subroutine WriteAtomicHamil(this, f)
    use AtLibrary, only: find
    class(AtomicHamil), intent(inout) :: this
    character(*), intent(in) :: f

    if(find(f,".kshell.snt")) then
      if(this%ms%sps%basis == "HO" .or. this%ms%sps%basis == "ho") then
        call write_hamil_kshell_snt(this, f)
        return
      end if
    end if

    if(find(f,".snt")) then
      call write_hamil_snt(this, f)
      return
    end if

    write(*,"(a)") "Unknown file format, writing file with snt format"
    call write_hamil_snt(this, f)

  end subroutine WriteAtomicHamil

  subroutine write_hamil_snt(this, f)
    use AtLibrary, only:hc, m_e
    class(AtomicHamil), intent(inout) :: this
    character(*), intent(in) :: f
    type(EleTwoBodySpace), pointer :: ms
    type(EleOrbits), pointer :: sps
    integer :: a, b, c, d, bra, ket, J, ch
    integer :: cnt, n
    integer :: wunit = 20
    ms => this%ms
    sps => ms%sps

    open(wunit, file=f, action="write")
#ifdef VERSION
    write(wunit, "(2a)") "# Generated by NuHamil (Tokyo code), ", trim(VERSION)
#endif
    write(wunit, "(3a)") "# Atomic Hamiltonian with ", trim(sps%basis), " basis"
    write(wunit, "(a,f6.4,a)") "# length scale is a_0 / ", ms%zeta, " with bohr radius a_0"
    write(wunit, "(a,f8.4,a)") "# Corresponding hw is ", &
        & hc **2 * this%ms%zeta**2 / (m_e*1.d3 * a_0**2), " eV"
    write(wunit, "(a)") "# Definition of model space"
    write(wunit, "(a)") "# numbers of proton orbit, neutron orbit, proton core, neutron core"
    write(wunit, '(4i5)') sps%norbs, 0, 0, 0
    write(wunit, "(a)") "# index, n,  l,  j,  z,  e"
    do bra = 1, sps%norbs
      write(wunit, '(6i5)') bra, sps%orb(bra)%n, sps%orb(bra)%l, sps%orb(bra)%j, -1, sps%orb(bra)%e
    end do

    cnt = 0
    do bra = 1, sps%norbs
      do ket = 1, bra
        if(sps%orb(bra)%l /= sps%orb(ket)%l) cycle
        if(sps%orb(bra)%j /= sps%orb(ket)%j) cycle
        cnt = cnt + 1
      end do
    end do

    write(wunit, '(a)') '####### one-body term'
    write(wunit, '(i5, i3)') cnt, 0
    write(wunit, '(a)') '### a b <a|t|b> <a|1/r|b> <a|rel_corr|b> <a|Darwin|b> <a|LS|b> (a.u.), <a|b>'
    do bra = 1, sps%norbs
      do ket = 1, bra
        if(sps%orb(bra)%l /= sps%orb(ket)%l) cycle
        if(sps%orb(bra)%j /= sps%orb(ket)%j) cycle
        write(wunit,'(2i5,6es20.10)') bra,ket, this%kinetic%m(bra,ket), this%potential%m(bra,ket), &
          & this%kinetic_p4%m(bra,ket), this%Darwin_term%m(bra,ket), &
          & this%LS_term%m(bra,ket), this%S%m(bra,ket)
        !write(wunit,'(2i5,3f18.8)') bra,ket, this%kinetic%m(bra,ket)*Eh, this%potential%m(bra,ket)*Eh, this%S%m(bra,ket)
      end do
    end do

    cnt = 0
    do ch = 1, ms%NChan
      cnt = cnt + ms%jp(ch)%n_state * (ms%jp(ch)%n_state + 1) / 2
    end do

    write(wunit, '(a)') '##### two-body term #####'
    write(wunit, '(i15, i3)') cnt, 0
    write(wunit, '(a)') '## a, b, c, d, J, coulomb, darwin, spin-contact, spin-orbit, orbit-orbit, spin-dipole (a.u.) '
    do ch = 1, ms%NChan
      j = ms%jp(ch)%j
      n = ms%jp(ch)%n_state
      if(n < 1) cycle
      do bra = 1, n
        a = ms%jp(ch)%n2label1(bra)
        b = ms%jp(ch)%n2label2(bra)
        do ket = 1, bra
          c = ms%jp(ch)%n2label1(ket)
          d = ms%jp(ch)%n2label2(ket)
          write(wunit, '(5i5, 6es20.10)') a, b, c, d, J, this%ee_coulomb(ch,ch)%m(bra, ket), &
              & this%ee_darwin(ch,ch)%m(bra,ket), this%ee_spin_contact(ch,ch)%m(bra,ket), &
              & this%ee_spin_orbit(ch,ch)%m(bra,ket), this%ee_orbit_orbit(ch,ch)%m(bra,ket), &
              & this%ee_spin_dipole(ch,ch)%m(bra,ket)
        end do
      end do
    end do

    close(wunit)
  end subroutine write_hamil_snt

  subroutine write_hamil_kshell_snt(this, f)
    !use SingleParticleState
    class(AtomicHamil), intent(inout) :: this
    character(*), intent(in) :: f
    type(EleTwoBodySpace), pointer :: ms
    type(EleOrbits), pointer :: sps
    type(EleOrbits) :: sps_k
    integer :: a, b, c, d, bra, ket, J, ch
    integer :: a_k, b_k, c_k, d_k
    integer :: cnt, n
    integer :: wunit = 20
    ms => this%ms
    sps => ms%sps

    call sps_k%init(sps%emax, mode="kshell")
    open(wunit, file=f, action="write")
    write(wunit, '(4i5)') sps_k%norbs/2, sps_k%norbs/2, 0, 0
    do bra = 1, sps_k%norbs
      write(wunit, '(5i5)') bra, sps_k%orb(bra)%n, sps_k%orb(bra)%l, sps_k%orb(bra)%j, -1
    end do

    cnt = 0
    do bra = 1, sps_k%norbs
      do ket = 1, bra
        if(abs(sps_k%orb(bra)%n - sps_k%orb(ket)%n) > 1) cycle
        if(sps_k%orb(bra)%l /= sps_k%orb(ket)%l) cycle
        if(sps_k%orb(bra)%j /= sps_k%orb(ket)%j) cycle
        cnt = cnt + 1
      end do
    end do

    write(wunit, '(a)') '####### one-body term'
    write(wunit, '(i5, i3)') cnt, 0
    write(wunit, '(a)') '### a, b, <a|t|b>'
    do bra = 1, sps_k%norbs
      do ket = 1, bra
        if(sps_k%orb(bra)%l /= sps_k%orb(ket)%l) cycle
        if(sps_k%orb(bra)%j /= sps_k%orb(ket)%j) cycle
        write(wunit,'(2i5,f18.8)') bra,ket, &
            & this%kinetic%m(bra,ket) - 1.d0 / dble(this%ms%zeta) * this%potential%m(bra,ket)
      end do
    end do

    cnt = 0
    do ch = 1, ms%NChan
      cnt = cnt + ms%jp(ch)%n_state * (ms%jp(ch)%n_state + 1) / 2
    end do

    write(wunit, '(a)') '##### two-body term #####'
    write(wunit, '(i15, i3)') cnt, 0
    write(wunit, '(a)') '## a, b, c, d, J, <ab|v|cd>'
    do ch = 1, ms%NChan
      j = ms%jp(ch)%j
      n = ms%jp(ch)%n_state
      if(n < 1) cycle
      do bra = 1, n
        a = ms%jp(ch)%n2label1(bra)
        b = ms%jp(ch)%n2label2(bra)
        a_k = sps_k%nlj2idx(sps%orb(a)%n, sps%orb(a)%l, sps%orb(a)%j)
        b_k = sps_k%nlj2idx(sps%orb(b)%n, sps%orb(b)%l, sps%orb(b)%j)
        do ket = 1, bra
          c = ms%jp(ch)%n2label1(ket)
          d = ms%jp(ch)%n2label2(ket)
          c_k = sps_k%nlj2idx(sps%orb(c)%n, sps%orb(c)%l, sps%orb(c)%j)
          d_k = sps_k%nlj2idx(sps%orb(d)%n, sps%orb(d)%l, sps%orb(d)%j)
          write(wunit, '(5i5, 1f16.8)') a_k, b_k, c_k, d_k, J, this%ee_coulomb(ch,ch)%m(bra, ket)
        end do
      end do
    end do

    close(wunit)
  end subroutine write_hamil_kshell_snt

  function kinetic_energy_laguerre_func(n1, n2, l) result(r)
    use AtLibrary, only: ln_gamma
    integer, intent(in) :: n1, n2, l
    real(8) :: ln_fact
    real(8) :: r

    r = 0.d0
    if(n1 > n2) then
      ln_fact = 0.5d0*ln_gamma(dble(n1+1)) + 0.5d0*ln_gamma(dble(n2+2*l+3)) - &
          &     0.5d0*ln_gamma(dble(n2+1)) - 0.5d0*ln_gamma(dble(n1+2*l+3))
      r = - dble(4*n2+4*l+6) / dble(2*l+3) * exp(ln_fact)
    end if

    if(n1 == n2) then
      r = - dble(4*n1+2*l+3) / dble(2*l+3)
    end if

    if(n1 < n2) then
      ln_fact = 0.5d0*ln_gamma(dble(n2+1)) + 0.5d0*ln_gamma(dble(n1+2*l+3)) - &
          &     0.5d0*ln_gamma(dble(n1+1)) - 0.5d0*ln_gamma(dble(n2+2*l+3))
      r = - dble(4*n1+4*l+6) / dble(2*l+3) * exp(ln_fact)
    end if
    r = -0.5d0 * r
  end function kinetic_energy_laguerre_func

  function coulomb_laguerre_func(n1, n2, l) result(r)
    use AtLibrary, only: ln_gamma
    integer, intent(in) :: n1, n2, l
    real(8) :: ln_fact
    real(8) :: r
    r = 0.d0
    if(n1 > n2) then
      ln_fact = 0.5d0*ln_gamma(dble(n1+1)) + 0.5d0*ln_gamma(dble(n2+2*l+3)) - &
          & 0.5d0*ln_gamma(dble(n2+1)) - 0.5d0*ln_gamma(dble(n1+2*l+3))
      r = exp(ln_fact) / dble(l+1)
    end if

    if(n1 <= n2) then
      ln_fact = 0.5d0*ln_gamma(dble(n2+1)) + 0.5d0*ln_gamma(dble(n1+2*l+3)) - &
          & 0.5d0*ln_gamma(dble(n1+1)) - 0.5d0*ln_gamma(dble(n2+2*l+3))
      r = exp(ln_fact) / dble(l+1)
    end if
  end function coulomb_laguerre_func
end module AtomicHamiltonian
