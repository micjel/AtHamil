module TwoBodyCoulomb
  use OneBodyTerms
  use ElectronTwoBodySpace
  implicit none

  private :: ee_interaction
  private :: finalize_fnl
  private :: initialize_fnl_laguerre
  private :: initialize_fnl_ho
  private :: initialize_fnl_hydrogen
  private :: get_fnl

  type, private :: FL
    real(8), allocatable :: F(:)
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
    procedure :: initialize_fnl_laguerre
    procedure :: initialize_fnl_hydrogen
    procedure :: initialize_fnl_ho
    procedure :: finalize_fnl
    procedure :: get_fnl
  end type Fnl
  type(Fnl), private :: Fintegral
  real(8), private, allocatable :: rmesh(:), rwmesh(:)
  real(8), private, allocatable :: rnl(:,:,:)
contains

  subroutine set_ee_coulomb_laguerre( this, ms, NMesh, rmax )
    use AtLibrary, only: gauss_legendre, laguerre_radial_wf_norm, &
        & fixed_point_quadrature, ln_gamma, laguerre
    type(AtomicHamilChan), intent(inout) :: this(:,:)
    type(EleTwoBodySpace), intent(in) :: ms
    integer, intent(in) :: NMesh
    integer :: n, l, i
    real(8) :: rmax, time
    write(*,*)
    write(*,"(a)") " Electron-electron Coulomb: "
#ifdef gauss_laguerre
    write(*,"(a)") "Using Gauss-Laguerre quadrature"
    call fixed_point_quadrature("laguerre", NMesh, rmesh, rwmesh, weight_renorm=.false., &
        & a_in=0.d0, b_in=2.d0, alpha_in=0.d0)
#else
    write(*,"(a)") "Using Gauss-Legendre quadrature"
    call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
#endif /* gauss_laguerre */
    allocate(rnl(NMesh, 0:ms%sps%emax, 0:ms%sps%lmax))
    rnl(:,:,:) = 0.d0
    do n = 0, ms%sps%emax
      do l = 0, ms%sps%lmax
        do i = 1, NMesh
#ifdef gauss_laguerre
          rnl(i,n,l) = exp( 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+2*l+3))) * &
              & laguerre(n,dble(2*l+2),2.d0*rmesh(i)) * (2.d0*rmesh(i))**(l+1) * sqrt(2.d0)
#else
          rnl(i,n,l) = laguerre_radial_wf_norm(n, l, 1.d0, rmesh(i))
#endif
        end do
      end do
    end do
    call norm_check(rnl, rwmesh, "lo")
    time = omp_get_wtime()
    call Fintegral%initialize_fnl_laguerre( ms )
    write(*,"(a,f12.6,a)") " Stored radial integral: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
    call set_ee_coulomb_term( this, ms )
    write(*,"(a,f12.6,a)") " Calculated ee Coulomb: ", omp_get_wtime() - time, " sec"
    call Fintegral%finalize_fnl()
  end subroutine set_ee_coulomb_laguerre

  subroutine set_ee_coulomb_ho( this, ms, NMesh, rmax )
    use AtLibrary, only: gauss_legendre, ho_radial_wf_norm
    type(AtomicHamilChan), intent(inout) :: this(:,:)
    type(EleTwoBodySpace), intent(in) :: ms
    integer, intent(in) :: NMesh
    integer :: n, l, i
    real(8) :: rmax, time

    write(*,*)
    write(*,"(a)") " Electron-electron Coulomb: "
    write(*,"(a)") "Using Gauss-Legendre quadrature"
    call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
    allocate(rnl(NMesh, 0:ms%sps%emax/2, 0:ms%sps%emax))
    rnl(:,:,:) = 0.d0
    do n = 0, ms%sps%emax/2
      do l = 0, ms%sps%emax
        do i = 1, NMesh
          rnl(i,n,l) = ho_radial_wf_norm(n, l, 1.d0, rmesh(i))
        end do
      end do
    end do
    call norm_check(rnl, rwmesh, "HO")
    time = omp_get_wtime()
    call Fintegral%initialize_fnl_ho( ms )
    write(*,"(a,f12.6,a)") " Stored radial integral: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
    call set_ee_coulomb_term( this, ms )
    write(*,"(a,f12.6,a)") " Calculated ee Coulomb: ", omp_get_wtime() - time, " sec"
    call Fintegral%finalize_fnl()
  end subroutine set_ee_coulomb_ho

  subroutine set_ee_coulomb_hydrogen( this, ms, NMesh, rmax )
    use AtLibrary, only: gauss_legendre, hydrogen_radial_wf_norm
    type(AtomicHamilChan), intent(inout) :: this(:,:)
    type(EleTwoBodySpace), intent(in) :: ms
    integer, intent(in) :: NMesh
    integer :: n, l, i
    real(8) :: rmax, time
    write(*,*)
    write(*,"(a)") " Electron-electron Coulomb: "
    write(*,"(a)") "Using Gauss-Legendre quadrature"
    call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
    allocate(rnl(NMesh, 1:ms%sps%emax, 0:ms%sps%lmax))
    rnl(:,:,:) = 0.d0
    do n = 1, ms%sps%emax
      do l = 0, min(n-1,ms%sps%lmax)
        do i = 1, NMesh
          rnl(i,n,l) = hydrogen_radial_wf_norm(n,l,1.d0,rmesh(i))
        end do
      end do
    end do
    call norm_check(rnl, rwmesh, "AO")
    time = omp_get_wtime()
    call Fintegral%initialize_fnl_hydrogen( ms )
    write(*,"(a,f12.6,a)") " Stored radial integral: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
    call set_ee_coulomb_term( this, ms )
    write(*,"(a,f12.6,a)") " Calculated ee Coulomb: ", omp_get_wtime() - time, " sec"
    call Fintegral%finalize_fnl()
  end subroutine set_ee_coulomb_hydrogen

  subroutine set_ee_coulomb_term( this, ms )
    type(AtomicHamilChan), intent(inout) :: this(:,:)
    type(EleTwoBodySpace), intent(in), target :: ms
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: bra, ket, a, b, c, d, J
    integer :: ch
    real(8) :: rabcd, rbacd, rabdc, rbadc, norm
    sps => ms%sps
    do ch = 1, ms%NChan
      J = ms%jp(ch)%j
      !$omp parallel
      !$omp do private(bra, a, b, oa, ob, ket, c, d, oc, od, norm, &
      !$omp &          rabcd, rbacd, rabdc, rbadc)
      do bra = 1, ms%jp(ch)%n_state
        a = ms%jp(ch)%n2label1(bra)
        b = ms%jp(ch)%n2label2(bra)
        oa => sps%orb(a)
        ob => sps%orb(b)
        do ket = 1, bra
          c = ms%jp(ch)%n2label1(ket)
          d = ms%jp(ch)%n2label2(ket)
          oc => sps%orb(c)
          od => sps%orb(d)

          norm = 1.d0
          if(a == b) norm = norm * dsqrt(0.5d0)
          if(c == d) norm = norm * dsqrt(0.5d0)
          rabcd = ee_interaction(oa, ob, oc, od, J)
          rbacd = ee_interaction(ob, oa, oc, od, J) * (-1.d0)**((oa%j+ob%j)/2-J-1)
          rabdc = ee_interaction(oa, ob, od, oc, J) * (-1.d0)**((oc%j+od%j)/2-J-1)
          rbadc = ee_interaction(ob, oa, od, oc, J) * (-1.d0)**((oa%j+ob%j+oc%j+od%j)/2)
          this(ch,ch)%m(bra,ket) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
          this(ch,ch)%m(ket,bra) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
  end subroutine set_ee_coulomb_term

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
      integral = Fintegral%get_fnl(oa%n,oa%l,ob%n,ob%l,oc%n,oc%l,od%n,od%l,L)
      r = r + integral * sjs(oa%j, ob%j, 2*J, od%j, oc%j, 2*L) * &
          & tjs(oa%j, 2*L, oc%j, -1, 0, 1) * &
          & tjs(ob%j, 2*L, od%j, -1, 0, 1)
    end do
    r = r * dsqrt(dble(oa%j+1) *dble(ob%j+1) * dble(oc%j+1) * dble(od%j+1)) * &
        &   (-1.d0) ** ((oa%j + oc%j) /2 + J)
  end function ee_interaction

  subroutine finalize_fnl(this)
    class(Fnl), intent(inout) :: this
    integer :: bra, ket
    do bra = 1, this%nidx
      do ket = 1, this%nidx
        deallocate(this%bk(bra,ket)%F)
      end do
    end do
    deallocate(this%bk)
    deallocate(this%idx)
    deallocate(this%n1)
    deallocate(this%n2)
    deallocate(this%l1)
    deallocate(this%l2)
  end subroutine finalize_fnl

  subroutine initialize_fnl_laguerre(this, ms)
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
    integer :: NMesh
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
        l1max = lmax
        l2max = lmax
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
        l1max = lmax
        l2max = lmax
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
        this%bk(bra,ket)%F(:) = 0.d0
      end do
    end do

    rmax_ = maxval(rmesh)
    NMesh = size(rmesh)
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
    deallocate(inter)
  end subroutine initialize_fnl_laguerre

  subroutine initialize_fnl_ho(this, ms)
    use AtLibrary, only: gauss_legendre, ho_radial_wf_norm
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
    integer :: NMesh
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
        l1max = lmax
        l2max = lmax
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
        l1max = lmax
        l2max = lmax
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
        this%bk(bra,ket)%F(:) = 0.d0
      end do
    end do

    rmax_ = maxval(rmesh)
    NMesh = size(rmesh)
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
          Rnl1(i1) = ho_radial_wf_norm(n1,l1,1.d0,r_1(i1)) * ho_radial_wf_norm(n3,l3,1.d0,r_1(i1)) * w_1(i1)
        end do
        do i1 = 1, n_2
          Rnl2(i1) = ho_radial_wf_norm(n1,l1,1.d0,r_2(i1)) * ho_radial_wf_norm(n3,l3,1.d0,r_2(i1)) * w_2(i1)
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
    deallocate(inter)
  end subroutine initialize_fnl_ho

  subroutine initialize_fnl_hydrogen(this, ms)
    use AtLibrary, only: gauss_legendre, hydrogen_radial_wf_norm
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
    integer :: NMesh
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
        l1max = min(n1 - 1, lmax)
        l2max = min(n2 - 1, lmax)
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
        l1max = min(n1 - 1, lmax)
        l2max = min(n2 - 1, lmax)
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
        this%bk(bra,ket)%F(:) = 0.d0
      end do
    end do

    rmax_ = maxval(rmesh)
    NMesh = size(rmesh)
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
          Rnl1(i1) = hydrogen_radial_wf_norm(n1,l1,1.d0,r_1(i1)) * hydrogen_radial_wf_norm(n3,l3,1.d0,r_1(i1)) * w_1(i1)
        end do
        do i1 = 1, n_2
          Rnl2(i1) = hydrogen_radial_wf_norm(n1,l1,1.d0,r_2(i1)) * hydrogen_radial_wf_norm(n3,l3,1.d0,r_2(i1)) * w_2(i1)
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
    deallocate(inter)
  end subroutine initialize_fnl_hydrogen

  function get_fnl(this, n1, l1, n2, l2, n3, l3, n4, l4, L) result(r)
    class(Fnl), intent(in) :: this
    integer, intent(in) :: n1, l1, n2, l2, n3, l3, n4, l4, L
    integer :: bra, ket
    real(8) :: r
    bra = this%idx(n1,l1,n3,l3)
    ket = this%idx(n2,l2,n4,l4)
    r = this%bk(bra,ket)%F(L)
  end function get_fnl
end module TwoBodyCoulomb
