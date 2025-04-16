module TwoBodySpinOrbit
  use OneBodyTerms
  use TwoBodyTerms
  use ElectronTwoBodySpace
  implicit none

  private :: laguerre_function_weight
  private :: d_laguerre_function_weight
  private :: set_ee_spin_orbit_term
  private :: ee_spin_orbit_interaction
  private :: C_function
  private :: CS_function
  private :: RS_function
  private :: R_function
  private :: R_function_l
  private :: finalize_fnl
  private :: initialize_fnl_laguerre
  private :: set_term_A
  private :: set_term_B
  private :: set_term_C
  private :: set_term_D
  private :: get_fnl
  private :: r_power_1_2
  private :: r_power_2_1

  type, private :: FL
    real(8), allocatable :: F(:,:)
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
    procedure :: finalize_fnl
    procedure :: get_fnl
  end type Fnl
  type(Fnl), private :: Fintegral
  real(8), private, allocatable :: rmesh(:), rwmesh(:)
  real(8), private, allocatable :: rnl(:,:,:)
contains

  subroutine set_ee_spin_orbit_laguerre( this, ms, NMesh, rmax )
    use AtLibrary, only: gauss_legendre, laguerre_radial_wf_norm, &
        & fixed_point_quadrature, ln_gamma, laguerre, d_laguerre_radial_wf_norm
    type(TwoBodyOperator), intent(inout) :: this
    type(EleTwoBodySpace), intent(in) :: ms
    integer, intent(in) :: NMesh
    integer :: n, l, i
    real(8) :: rmax, time
    write(*,*)
    write(*,"(a)") " Electron-electron spin-orbit: "
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
          rnl(i,n,l) = laguerre_radial_wf_norm_glmesh(n, dble(l), 1.d0, rmesh(i))
#else
          rnl(i,n,l) = laguerre_radial_wf_norm(n, dble(l), 1.d0, rmesh(i))
#endif
        end do
      end do
    end do
    call norm_check(rnl, rwmesh, "lo")
    time = omp_get_wtime()
    call Fintegral%initialize_fnl_laguerre( ms )
    write(*,"(a,f12.6,a)") " Stored radial integral: ", omp_get_wtime() - time, " sec"
    time = omp_get_wtime()
    call set_ee_spin_orbit_term( this, ms )
    write(*,"(a,f12.6,a)") " Calculated ee spin-orbit: ", omp_get_wtime() - time, " sec"
    call Fintegral%finalize_fnl()
  end subroutine set_ee_spin_orbit_laguerre

  function laguerre_function_weight(n, l, r) result(s)
    use AtLibrary, only: ln_gamma, laguerre
    integer, intent(in) :: n
    real(8), intent(in) :: l, r
    real(8) :: s
    s = exp( 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+3)+2.d0*l) + (l+1.d0)*log(2.d0*r) ) * &
        & laguerre(n,2.d0*l+2.d0,2.d0*r) * sqrt(2.d0)
  end function laguerre_function_weight

  function d_laguerre_function_weight(n, l, zeta, r) result(s)
    integer, intent(in) :: n
    real(8), intent(in) :: l, zeta, r
    real(8) :: s
    real(8) :: x
    x = 2.d0 * r / zeta
    s = ( (0.5d0 * zeta * l) / x - 0.25d0 * zeta ) * laguerre_function_weight(n,l,r)
    if(n==0) return
    s = s - 0.5d0 * zeta * sqrt(dble(n)) / sqrt(x) * laguerre_function_weight(n-1,l+0.5d0,r)
  end function d_laguerre_function_weight

  subroutine set_ee_spin_orbit_term( this, ms )
    type(TwoBodyOperator), intent(inout) :: this
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
          rabcd = ee_spin_orbit_interaction(oa, ob, oc, od, J)
          rbacd = ee_spin_orbit_interaction(ob, oa, oc, od, J) * (-1.d0)**((oa%j+ob%j)/2-J-1)
          rabdc = ee_spin_orbit_interaction(oa, ob, od, oc, J) * (-1.d0)**((oc%j+od%j)/2-J-1)
          rbadc = ee_spin_orbit_interaction(ob, oa, od, oc, J) * (-1.d0)**((oa%j+ob%j+oc%j+od%j)/2)
          this%MatCh(ch,ch)%m(bra,ket) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
          this%MatCh(ch,ch)%m(ket,bra) = 0.5d0 * norm * (rabcd + rbacd + rabdc + rbadc)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end do
  end subroutine set_ee_spin_orbit_term

  function ee_spin_orbit_interaction(oa, ob, oc, od, J) result(r)
    use AtLibrary, only: tjs, sjs, alpha, triag
    type(EleSingleParticleOrbit), intent(in) :: oa, ob, oc, od
    integer, intent(in) :: J
    real(8) :: r
    r = ee_spin_orbit_interaction_(oa,ob,oc,od,J) - 0.5d0* &
        & ee_spin_orbit_interaction_(ob,oa,od,oc,J)
  end function ee_spin_orbit_interaction

  function ee_spin_orbit_interaction_(oa, ob, oc, od, J) result(r)
    use AtLibrary, only: tjs, sjs, alpha, triag
    type(EleSingleParticleOrbit), intent(in) :: oa, ob, oc, od
    integer, intent(in) :: J
    real(8) :: r
    integer :: Lmin, Lmax, L
    real(8) :: I1, I2, I3, I4
    Lmin = max(abs(oa%j-oc%j), abs(ob%j-od%j))/2
    Lmax = min(   (oa%j+oc%j),    (ob%j+od%j))/2
    r = 0.d0
    do L = Lmin, Lmax
      if( triag( oa%l, oc%l, L) ) cycle
      if( triag( ob%l, od%l, L) ) cycle
      if( mod(oa%l+oc%l+L,2) == 1) cycle
      if( mod(ob%l+od%l+L,2) == 1) cycle
      I1 = Fintegral%get_fnl( oa%n, oa%l, ob%n, ob%l, oc%n, oc%l, od%n, od%l, L, 1 )
      I2 = Fintegral%get_fnl( oa%n, oa%l, ob%n, ob%l, oc%n, oc%l, od%n, od%l, L, 2 )
      I3 = Fintegral%get_fnl( oa%n, oa%l, ob%n, ob%l, oc%n, oc%l, od%n, od%l, L, 3 )
      I4 = Fintegral%get_fnl( oa%n, oa%l, ob%n, ob%l, oc%n, oc%l, od%n, od%l, L, 4 )

      r = r - sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & I1 * sqrt( dble(L*(L+1))/dble(2*L+1) ) * &
          & CS_function(L,L,oa%l,oa%j,oc%l,oc%j) * C_function(L,ob%l,ob%j,od%l,od%j)
      r = r + sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) / ( dble(2*L+1)) * &
          & I2 * RS_function(L,L,L,oa%l,oa%j,oc%l,oc%j) * C_function(L,ob%l,ob%j,od%l,od%j)
      r = r + sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & dble(L+2) / sqrt( dble( (2*L+1) * (2*L+3) )) * &
          & I3 * RS_function(L,L+1,L,oa%l,oa%j,oc%l,oc%j) * C_function(L,ob%l,ob%j,od%l,od%j)
      r = r + sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & sqrt( dble((L+1)*(L+2))) / sqrt( dble( (2*L+1) * (2*L+3) )) * &
          & I3 * RS_function(L+2,L+1,L,oa%l,oa%j,oc%l,oc%j) * C_function(L,ob%l,ob%j,od%l,od%j)
      r = r - 0.5d0 * sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & sqrt( dble(L*(L+1))/dble(2*L+1) ) * &
          & I1 * C_function(L,oa%l,oa%j,oc%l,oc%j) * CS_function(L,L,ob%l,ob%j,od%l,od%j)
      r = r + 0.5d0 * sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & 1.d0 / ( dble(2*L+1)) * &
          & I2 * R_function(L,L,oa%l,oa%j,oc%l,oc%j) * CS_function(L,L,ob%l,ob%j,od%l,od%j)
      r = r + 0.5d0 * sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & dble(L+2) / sqrt( dble( (2*L+1) * (2*L+3) )) * &
          & I3 * R_function(L,L+1,oa%l,oa%j,oc%l,oc%j) * CS_function(L,L+1,ob%l,ob%j,od%l,od%j)
      r = r + 0.5d0 * sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & sqrt( dble((L+1)*(L+2))) / sqrt( dble( (2*L+1) * (2*L+3) )) * &
          & I3 * R_function(L+2,L+1,oa%l,oa%j,oc%l,oc%j) * CS_function(L,L+1,ob%l,ob%j,od%l,od%j)
      if( L == 0 ) cycle
      r = r - sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          &  dble((L-1)) / sqrt( dble( (2*L+1) * (2*L-1) )) * &
          & I4 * RS_function(L-1,L,L,oa%l,oa%j,oc%l,oc%j) * C_function(L,ob%l,ob%j,od%l,od%j)
      r = r - 0.5d0 * sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & dble((L-1))/sqrt( dble( (2*L+1)*(2*L+3) )) * &
          & I4 * R_function(L-1,L,oa%l,oa%j,oc%l,oc%j) * CS_function(L,L-1,ob%l,ob%j,od%l,od%j)
      if( L == 1 ) cycle
      r = r - sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & sqrt( dble((L-1)*L)) / sqrt( dble( (2*L+1) * (2*L-1) )) * &
          & I4 * RS_function(L-2,L-1,L,oa%l,oa%j,oc%l,oc%j) * C_function(L,ob%l,ob%j,od%l,od%j)
      r = r - 0.5d0 * sjs( oa%j, ob%j, 2*J, od%j, oc%j, 2*L ) * &
          & sqrt( dble((L-1)*L)) / sqrt( dble( (2*L+1) * (2*L+3) )) * &
          & I4 * R_function(L-2,L-1,oa%l,oa%j,oc%l,oc%j) * CS_function(L,L-1,ob%l,ob%j,od%l,od%j)
    end do
    r = r * (-1.d0)**((ob%j + oc%j) /2 + J) * 4*pi / alpha**2
  end function ee_spin_orbit_interaction_

  function C_function(L, la, ja, lc, jc) result(r)
    ! ( a || Y^L || c )
    use AtLibrary, only: tjs, triag
    integer, intent(in) :: L, la, ja, lc, jc
    real(8) :: r
    r = 0.d0
    if( mod(la+lc+L,2) == 1) return
    r = (-1.d0)**( (ja-1)/2 ) * sqrt( dble( (ja+1)*(jc+1)*(2*L+1) )) * &
        & tjs( ja, 2*L, jc, -1, 0, 1 )
  end function C_function

  function CS_function(K, L, la, ja, lc, jc) result(r)
    ! ( a || [ Y^K S ]^L || c )
    use AtLibrary, only: tjs, snj
    integer, intent(in) :: K, L, la, ja, lc, jc
    real(8) :: r
    r = 0.d0
    if( mod(la+lc+K,2) == 1) return
    r = (-1)**la * sqrt( 1.5d0*dble(ja+1)*dble(2*L+1)*dble(jc+1)*dble(2*la+1)*dble(2*lc+1)*dble(2*K+1) ) * &
        & snj( 2*la, 2*lc, 2*K, 1, 1, 2, ja, jc, 2*L ) * &
        & tjs( 2*la, 2*K, 2*lc, 0, 0, 0)
  end function CS_function

  function RS_function(G, K, L, la, ja, lc, jc) result(r)
    use AtLibrary, only: snj
    integer, intent(in) :: G, K, L, la, ja, lc, jc
    real(8) :: r
    r = 0.d0
    if( mod(la+lc+G,2) == 1 ) return
    if( lc == 0 ) return
    r = sqrt( dble(ja+1)*dble(2*L+1)*dble(jc+1)*1.5d0 ) * &
        & snj(2*la, 2*lc, 2*K, 1, 1, 2, ja, jc, 2*L) * &
        & R_function_l(G, K, la, lc)
  end function RS_function

  function R_function_l(K, L, la, lc) result(r)
    ! ( la || [Y^K L]^L || lc )
    use AtLibrary, only: tjs, sjs, triag
    integer, intent(in) :: K, L, la, lc
    real(8) :: r
    r = 0.d0
    if( mod(la+lc+K,2) == 1) return
    if( lc==0 ) return
    r = (-1.d0)**(L+lc) * sqrt(dble(2*L+1)) * &
        & sjs(2*K, 2*L, 2, 2*lc, 2*lc, 2*la) * &
        & sqrt( dble( lc*(lc+1)*(2*lc+1) )) * &
        & sqrt( dble( (2*la+1)*(2*lc+1)*(2*K+1) )) * tjs(2*la, 2*K, 2*lc, 0, 0, 0)
  end function R_function_l

  function R_function(K, L, la, ja, lc, jc) result(r)
    ! ( a || [Y^K L]^L || c )
    use AtLibrary, only: tjs, sjs, triag
    integer, intent(in) :: K, L, la, ja, lc, jc
    real(8) :: r
    r = 0.d0
    if( mod(la+lc+K,2) == 1) return
    if( lc==0 ) return
    r = (-1.d0)**( la+lc+(jc+1)/2 ) * &
        & sqrt( dble( (ja+1)*(jc+1)*(2*la+1)*(2*lc+1)*(2*L+1)*(2*K+1) )) * &
        & sqrt( dble( lc * (lc+1) * (2*lc+1) )) * &
        & sjs( 2*la, ja, 1, jc, 2*lc, 2*L ) * &
        & sjs( 2*K, 2*L, 2, 2*lc, 2*lc, 2*la ) * &
        & tjs( 2*la, 2*K, 2*lc, 0, 0, 0 )
  end function R_function

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
    class(Fnl), intent(inout) :: this
    type(EleTwoBodySpace), intent(in) :: ms
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa
    integer :: nmin, nmax, lmin, lmax
    integer :: l1max, l2max
    integer :: n1, l1, n2, l2
    integer :: n3, l3, n4, l4
    integer :: bra, ket
    integer :: cnt, a
    real(8) :: zeta

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

        lmin = max(max(abs(l1-l3), abs(l2-l4))-1, 0)
        lmax = min(abs(l1+l3), abs(l2+l4))+1
        allocate(this%bk(bra,ket)%F(lmin:lmax,4))
        this%bk(bra,ket)%F(:,:) = 0.d0
      end do
    end do

    call set_term_A(this, zeta, lmax)
    call set_term_B(this, zeta, lmax)
    call set_term_C(this, zeta, lmax)
    call set_term_D(this, zeta, lmax)
  end subroutine initialize_fnl_laguerre

  subroutine set_term_A(this, zeta, llmax)
    use AtLibrary, only: gauss_legendre, laguerre_radial_wf_norm, d_laguerre_radial_wf_norm, triag
    type(Fnl), intent(inout) :: this
    integer, intent(in) :: llmax
    real(8), intent(in) :: zeta
    real(8) :: integral, r1, r2
    real(8), allocatable :: inter(:,:,:)
    real(8), allocatable, save :: r_1(:), w_1(:)
    real(8), allocatable, save :: r_2(:), w_2(:)
    real(8), allocatable, save :: dRnl1(:), dRnl2(:)
    real(8) :: rmax_
    integer :: NMesh
    integer :: n_1, n_2
    integer :: i1, i2, L, lmin, lmax
    integer :: n1, l1, n2, l2
    integer :: n3, l3, n4, l4
    integer :: bra, ket
    !$omp threadprivate(r_1, w_1, r_2, w_2, dRnl1, dRnl2)
    rmax_ = maxval(rmesh)
    NMesh = size(rmesh)
    allocate(inter(NMesh,0:2*llmax+1,this%nidx))
    inter(:,:,:) = 0.d0
    !$omp parallel
    !$omp do private(i2, r2, n_1, n_2, ket, n1, l1, n3, l3, i1, L, integral, r1)
    do i2 = 1, NMesh
      r2 = rmesh(i2) / zeta
      n_1 = max(int(rmesh(i2) / rmax_ * dble(NMesh)),NMesh/20)
      n_2 = NMesh - n_1
      call gauss_legendre(0.d0, rmesh(i2), r_1, w_1, n_1)
      call gauss_legendre(rmesh(i2), rmax_, r_2, w_2, n_2)
      allocate( DRnl1(n_1) )
      allocate( DRnl2(n_2) )
      do ket = 1, this%nidx
        n1 = this%n1(ket)
        l1 = this%l1(ket)
        n3 = this%n2(ket)
        l3 = this%l2(ket)
        do i1 = 1, n_1
          dRnl1(i1) = laguerre_radial_wf_norm(n1,dble(l1),1.d0,r_1(i1)) * &
              & d_laguerre_radial_wf_norm(n3,dble(l3),1.d0,r_1(i1)) * w_1(i1)
        end do
        do i1 = 1, n_2
          dRnl2(i1) = laguerre_radial_wf_norm(n1,dble(l1),1.d0,r_2(i1)) * &
              & d_laguerre_radial_wf_norm(n3,dble(l3),1.d0,r_2(i1)) * w_2(i1)
        end do

        do L = abs(l1-l3)-1, l1+l3+1
          if(L < 0) cycle
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle
          integral = 0.d0
          do i1 = 1, n_1 ! r1 < r2
            r1 = r_1(i1) / zeta
            integral = integral + r_power_1_2(r1,r2,L-1,L+1)*dRnl1(i1)
          end do
          do i1 = 1, n_2 ! r1 > r2
            r1 = r_2(i1) / zeta
            integral = integral + r_power_2_1(r1,r2,L,L+2)*dRnl2(i1)
          end do
          inter(i2,L,ket) = integral

        end do
      end do
      deallocate(r_1, w_1, r_2, w_2, DRnl1, DRnl2)
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
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle
          if( mod(l2+l4+L,2) == 1) cycle
          if( triag(l2,l4,L)) cycle
          integral = 0.d0
          do i2 = 1, NMesh
            integral = integral + inter(i2,L,ket) * rwmesh(i2) * &
                & rnl(i2,n1,l1) * rnl(i2,n3,l3)
          end do
          this%bk(bra,ket)%F(L,1) = integral
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(inter)
  end subroutine set_term_A

  subroutine set_term_B(this, zeta, llmax)
    use AtLibrary, only: gauss_legendre, laguerre_radial_wf_norm, d_laguerre_radial_wf_norm, triag
    type(Fnl), intent(inout) :: this
    integer, intent(in) :: llmax
    real(8), intent(in) :: zeta
    real(8) :: integral, r1, r2
    real(8), allocatable :: inter(:,:,:)
    real(8), allocatable, save :: r_1(:), w_1(:)
    real(8), allocatable, save :: r_2(:), w_2(:)
    real(8), allocatable, save :: Rnl1(:), Rnl2(:)
    real(8) :: rmax_
    integer :: NMesh
    integer :: n_1, n_2
    integer :: i1, i2, L, lmin, lmax
    integer :: n1, l1, n2, l2
    integer :: n3, l3, n4, l4
    integer :: bra, ket
    !$omp threadprivate(r_1, w_1, r_2, w_2, Rnl1, Rnl2)
    rmax_ = maxval(rmesh)
    NMesh = size(rmesh)
    allocate(inter(NMesh,0:2*llmax+1,this%nidx))
    inter(:,:,:) = 0.d0
    !$omp parallel
    !$omp do private(i2, r2, n_1, n_2, ket, n1, l1, n3, l3, i1, L, integral, r1)
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
          Rnl1(i1) =  laguerre_radial_wf_norm(n1,dble(l1),1.d0,r_1(i1)) * &
              & laguerre_radial_wf_norm(n3,dble(l3),1.d0,r_1(i1)) * w_1(i1)
        end do
        do i1 = 1, n_2
          Rnl2(i1) =  laguerre_radial_wf_norm(n1,dble(l1),1.d0,r_2(i1)) * &
              & laguerre_radial_wf_norm(n3,dble(l3),1.d0,r_2(i1)) * w_2(i1)
        end do

        do L = abs(l1-l3)-1, l1+l3+1
          if(L < 0) cycle
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle

          integral = 0.d0
          do i1 = 1, n_1 ! r1 > r2
            r1 = r_1(i1) / zeta
            integral = integral - dble(L) * r_power_2_1(r1,r2,L,L+3) * Rnl1(i1)
          end do
          do i1 = 1, n_2 ! r1 < r2
            r1 = r_2(i1) / zeta
            integral = integral + dble(L+1) * r_power_1_2(r1,r2,L-1,L+2) * Rnl2(i1)
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
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle
          if( mod(l2+l4+L,2) == 1) cycle
          if( triag(l2,l4,L)) cycle
          integral = 0.d0
          do i2 = 1, NMesh
            integral = integral + inter(i2,L,ket) * rwmesh(i2) * &
                & rnl(i2,n1,l1) * rnl(i2,n3,l3)
          end do
          this%bk(bra,ket)%F(L,2) = integral
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(inter)
  end subroutine set_term_B

  subroutine set_term_C(this, zeta, llmax)
    use AtLibrary, only: gauss_legendre, laguerre_radial_wf_norm, triag
    type(Fnl), intent(inout) :: this
    integer, intent(in) :: llmax
    real(8), intent(in) :: zeta
    real(8) :: integral, r1, r2
    real(8), allocatable :: inter(:,:,:)
    real(8), allocatable, save :: r_1(:), w_1(:)
    real(8), allocatable, save :: r_2(:), w_2(:)
    real(8), allocatable, save :: Rnl1(:), Rnl2(:)
    real(8) :: rmax_
    integer :: NMesh
    integer :: n_1, n_2
    integer :: i1, i2, L, lmin, lmax
    integer :: n1, l1, n2, l2
    integer :: n3, l3, n4, l4
    integer :: bra, ket
    !$omp threadprivate(r_1, w_1, r_2, w_2, Rnl1, Rnl2)
    rmax_ = maxval(rmesh)
    NMesh = size(rmesh)
    allocate(inter(NMesh,0:2*llmax+1,this%nidx))
    inter(:,:,:) = 0.d0
    !$omp parallel
    !$omp do private(i2, r2, n_1, n_2, ket, n1, l1, n3, l3, i1, L, integral, r1)
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
          Rnl1(i1) =  laguerre_radial_wf_norm(n1,dble(l1),1.d0,r_1(i1)) * &
              & laguerre_radial_wf_norm(n3,dble(l3),1.d0,r_1(i1)) * w_1(i1)
        end do
        do i1 = 1, n_2
          Rnl2(i1) =  0.d0
        end do

        do L = abs(l1-l3)-1, l1+l3+1
          if(L < 0) cycle
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle

          integral = 0.d0
          do i1 = 1, n_1 ! r1 > r2
            r1 = r_1(i1) / zeta
            integral = integral - r_power_2_1(r1,r2,L,L+3) * Rnl1(i1)
          end do
          do i1 = 1, n_2 ! r1 < r2
            r1 = r_2(i1) / zeta
            integral = integral - r_power_2_1(r1,r2,L,L+3) * Rnl2(i1)
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
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle
          if( mod(l2+l4+L,2) == 1) cycle
          if( triag(l2,l4,L)) cycle
          integral = 0.d0
          do i2 = 1, NMesh
            integral = integral + inter(i2,L,ket) * rwmesh(i2) * &
                & rnl(i2,n1,l1) * rnl(i2,n3,l3)
          end do
          this%bk(bra,ket)%F(L,3) = integral
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(inter)
  end subroutine set_term_C

  subroutine set_term_D(this, zeta, llmax)
    use AtLibrary, only: gauss_legendre, laguerre_radial_wf_norm, triag
    type(Fnl), intent(inout) :: this
    integer, intent(in) :: llmax
    real(8), intent(in) :: zeta
    real(8) :: integral, r1, r2
    real(8), allocatable :: inter(:,:,:)
    real(8), allocatable, save :: r_1(:), w_1(:)
    real(8), allocatable, save :: r_2(:), w_2(:)
    real(8), allocatable, save :: Rnl1(:), Rnl2(:)
    real(8) :: rmax_
    integer :: NMesh
    integer :: n_1, n_2
    integer :: i1, i2, L, lmin, lmax
    integer :: n1, l1, n2, l2
    integer :: n3, l3, n4, l4
    integer :: bra, ket
    !$omp threadprivate(r_1, w_1, r_2, w_2, Rnl1, Rnl2)
    rmax_ = maxval(rmesh)
    NMesh = size(rmesh)
    allocate(inter(NMesh,0:2*llmax+1,this%nidx))
    inter(:,:,:) = 0.d0
    !$omp parallel
    !$omp do private(i2, r2, n_1, n_2, ket, n1, l1, n3, l3, i1, L, integral, r1)
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
          Rnl1(i1) =  0.d0
        end do
        do i1 = 1, n_2
          Rnl2(i1) =  laguerre_radial_wf_norm(n1,dble(l1),1.d0,r_2(i1)) * &
              & laguerre_radial_wf_norm(n3,dble(l3),1.d0,r_2(i1)) * w_2(i1)
        end do

        do L = abs(l1-l3)-1, l1+l3+1
          if(L < 0) cycle
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle

          integral = 0.d0
          do i1 = 1, n_1 ! r1 > r2
            r1 = r_1(i1) / zeta
            integral = integral + r_power_1_2(r1,r2,L-1,L+2) * Rnl1(i1)
          end do
          do i1 = 1, n_2 ! r1 < r2
            r1 = r_2(i1) / zeta
            integral = integral + r_power_1_2(r1,r2,L-1,L+2) * Rnl2(i1)
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
          if( mod(l1+l3+L,2) == 1) cycle
          if( triag(l1,l3,L)) cycle
          if( mod(l2+l4+L,2) == 1) cycle
          if( triag(l2,l4,L)) cycle
          integral = 0.d0
          do i2 = 1, NMesh
            integral = integral + inter(i2,L,ket) * rwmesh(i2) * &
                & rnl(i2,n1,l1) * rnl(i2,n3,l3)
          end do
          this%bk(bra,ket)%F(L,4) = integral
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
    deallocate(inter)
  end subroutine set_term_D

  function get_fnl(this, n1, l1, n2, l2, n3, l3, n4, l4, L, i) result(r)
    class(Fnl), intent(in) :: this
    integer, intent(in) :: n1, l1, n2, l2, n3, l3, n4, l4, L, i
    integer :: bra, ket
    real(8) :: r
    bra = this%idx(n1,l1,n3,l3)
    ket = this%idx(n2,l2,n4,l4)
    r = this%bk(bra,ket)%F(L,i)
  end function get_fnl

  function r_power_1_2(r1,r2,l_num,l_den) result(r)
    real(8), intent(in) :: r1, r2
    integer, intent(in) :: l_num, l_den
    real(8) :: r
    r = exp( l_num*log( r1 ) - l_den*log( r2 ) )
  end function r_power_1_2

  function r_power_2_1(r1,r2,l_num,l_den) result(r)
    real(8), intent(in) :: r1, r2
    integer, intent(in) :: l_num, l_den
    real(8) :: r
    r = exp( l_num*log( r2 ) - l_den*log( r1 ) )
  end function r_power_2_1
end module TwoBodySpinOrbit
