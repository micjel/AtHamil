module OneBodyTerms
  use omp_lib
  use LinAlgLib
  use EleSingleParticleState

  implicit none
  type, extends(DMat) :: OneBodyOperator
    type(EleOrbits), pointer :: sps
    integer :: jrank = 0, prank = 1
    logical :: Zero = .true.
  contains
    procedure :: InitOneBodyOperator
    procedure :: FinOneBodyOperator
    procedure :: SetOneBodyOperator
    generic :: init => InitOneBodyOperator
    generic :: release => FinOneBodyOperator
    generic :: set => SetOneBodyOperator
  end type OneBodyOperator

  real(8), private, allocatable :: rmesh(:), rwmesh(:)
  real(8), private, allocatable :: pmesh(:), pwmesh(:)
  real(8), private, allocatable :: rnl(:,:,:)
  real(8), private, allocatable :: rnl_mom(:,:,:)
contains

  subroutine InitOneBodyOperator( this, sps, jrank, prank )
    class(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    integer, intent(in), optional :: jrank, prank
    this%sps => sps
    if(present( jrank )) this%jrank = jrank
    if(present( prank )) this%prank = prank
    call this%zeros( sps%norbs, sps%norbs )
  end subroutine InitOneBodyOperator

  subroutine FinOneBodyOperator( this )
    class(OneBodyOperator), intent(inout) :: this
    call this%fin()
    this%sps => null()
    this%jrank = 0
    this%prank = 1
  end subroutine FinOneBodyOperator

  subroutine SetOneBodyOperator( this, opname, zeta, NMesh, xmax )
    class(OneBodyOperator), intent(inout) :: this
    character(*), intent(in) :: opname
    real(8), intent(in) :: zeta
    real(8), intent(in), optional :: xmax
    integer, intent(in), optional :: NMesh
    real(8) :: xxmax = 0.d0
    integer :: NNMesh = 0

    if(present(xmax)) xxmax = xmax
    if(present(NMesh)) NNMesh = NMesh
    select case(opname)
    case("kinetic_laguerre")
      call set_kinetic_laguerre(this, this%sps, zeta)
    case("kinetic_ho")
      call set_kinetic_ho(this, this%sps, zeta)
    case("kinetic_hydrogen")
      if(.not. present( xmax ) ) stop "Input pmax needed"
      if(.not. present( NMesh ) ) stop "Input NMesh needed"
      call set_kinetic_hydrogen(this, this%sps, zeta, NNMesh, xxmax)
    case("coulomb_laguerre")
      if(.not. present( xmax ) ) stop "Input pmax needed"
      if(.not. present( NMesh ) ) stop "Input NMesh needed"
      call set_coulomb_laguerre(this, this%sps, zeta, NNMesh, xxmax)
    case("coulomb_ho")
      if(.not. present( xmax ) ) stop "Input pmax needed"
      if(.not. present( NMesh ) ) stop "Input NMesh needed"
      call set_coulomb_ho(this, this%sps, zeta, NNMesh, xxmax)
    case("coulomb_hydrogen")
      if(.not. present( xmax ) ) stop "Input pmax needed"
      if(.not. present( NMesh ) ) stop "Input NMesh needed"
      call set_coulomb_hydrogen(this, this%sps, zeta, NNMesh, xxmax)
    case("kinetic_correction_laguerre")
      if(.not. present( xmax ) ) stop "Input pmax needed"
      if(.not. present( NMesh ) ) stop "Input NMesh needed"
      call set_kinetic_correction_laguerre(this, this%sps, zeta, NNMesh, xxmax)
    case("Darwin_laguerre")
      call set_Darwin_laguerre(this, this%sps, zeta)
    case("spin_orbit_laguerre")
      if(.not. present( xmax ) ) stop "Input pmax needed"
      if(.not. present( NMesh ) ) stop "Input NMesh needed"
      call set_spin_orbit_laguerre(this, this%sps, zeta, NNMesh, xxmax)
    end select
  end subroutine SetOneBodyOperator

  subroutine set_kinetic_laguerre(this, sps, zeta)
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    real(8) :: r

    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        r = kinetic_energy_laguerre_func(oa%n,ob%n,oa%l) * zeta**2
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
  end subroutine set_kinetic_laguerre

  subroutine set_kinetic_ho(this, sps, zeta)
    use AtLibrary, only: hc, m_e
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    real(8) :: r, hw

    hw = hc **2 * zeta**2 / m_e*1.d3 ! in a.u.
    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        if(abs(oa%n - ob%n) > 1) cycle
        r = 0.d0
        if(oa%n == ob%n) r = dble(2 * oa%n + oa%l) + 1.5d0
        if(oa%n == ob%n-1) r = dsqrt(dble(oa%n + 1) * (dble(oa%n + oa%l) + 1.5d0))
        if(oa%n == ob%n+1) r = dsqrt(dble(ob%n + 1) * (dble(ob%n + ob%l) + 1.5d0))
        r = r * hw * 0.5d0
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
  end subroutine set_kinetic_ho

  subroutine set_kinetic_hydrogen(this, sps, zeta, NMesh, pmax)
    use AtLibrary, only: hc, m_e, gauss_legendre, hydrogen_radial_wf_mom_norm
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta, pmax
    integer, intent(in) :: NMesh
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    integer :: n, l, i
    real(8) :: r

    call gauss_legendre(0.d0, pmax, pmesh, pwmesh, NMesh)
    allocate(rnl_mom(NMesh, 1:sps%emax, 0:sps%lmax))
    rnl_mom(:,:,:) = 0.d0
    do n = 1, sps%emax
      do l = 0, min(n-1,sps%lmax)
        do i = 1, NMesh
          rnl_mom(i,n,l) = hydrogen_radial_wf_mom_norm(n,l,1.d0,pmesh(i))
        end do
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
          r = r + rnl_mom(i,oa%n,oa%l) * rnl_mom(i,ob%n,ob%l) * &
              & pwmesh(i) * (pmesh(i)*zeta)**2 * 0.5d0
        end do
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
    deallocate(rnl_mom)
    deallocate(pmesh)
    deallocate(pwmesh)
  end subroutine set_kinetic_hydrogen

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

  subroutine set_coulomb_laguerre(this, sps, zeta, NMesh, rmax)
    use AtLibrary, only: hc, m_e, gauss_legendre, &
        & fixed_point_quadrature, ln_gamma, laguerre, &
        & laguerre_radial_wf_norm
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta, rmax
    integer, intent(in) :: NMesh
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    integer :: n, l, i
    real(8) :: r

#ifdef gauss_laguerre
    call fixed_point_quadrature("laguerre", NMesh, rmesh, rwmesh, weight_renorm=.false., &
        & a_in=0.d0, b_in=2.d0, alpha_in=0.d0)
#else
    call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
#endif /* gauss_laguerre */
    allocate(rnl(NMesh, 0:sps%emax, 0:sps%lmax))
    rnl(:,:,:) = 0.d0
    do n = 0, sps%emax
      do l = 0, sps%lmax
        do i = 1, NMesh
#ifdef gauss_laguerre
          rnl(i,n,l) = exp( 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+2*l+3))) * &
              & laguerre(n,dble(2*l+2),2.d0*rmesh(i)) * (2.d0*rmesh(i))**(l+1) * sqrt(2.d0)
#else
          rnl(i,n,l) = laguerre_radial_wf_norm(n, l, 1.d0, rmesh(i))
#endif /* gauss_laguerre */
        end do
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
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
    deallocate(rnl)
    deallocate(rmesh, rwmesh)
  end subroutine set_coulomb_laguerre

  subroutine set_coulomb_HO(this, sps, zeta, NMesh, rmax)
    use AtLibrary, only: hc, m_e, gauss_legendre, &
        & ho_radial_wf_norm
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta, rmax
    integer, intent(in) :: NMesh
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    integer :: n, l, i
    real(8) :: r

    call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
    allocate(rnl(NMesh, 0:sps%emax/2, 0:sps%lmax))
    rnl(:,:,:) = 0.d0
    do n = 0, sps%emax/2
      do l = 0, sps%lmax
        do i = 1, NMesh
          rnl(i,n,l) = ho_radial_wf_norm(n, l, 1.d0, rmesh(i))
        end do
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
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
    deallocate(rnl)
    deallocate(rmesh, rwmesh)
  end subroutine set_coulomb_HO

  subroutine set_coulomb_hydrogen(this, sps, zeta, NMesh, rmax)
    use AtLibrary, only: hc, m_e, gauss_legendre, &
        & hydrogen_radial_wf_norm
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta, rmax
    integer, intent(in) :: NMesh
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    integer :: n, l, i
    real(8) :: r

    call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
    allocate(rnl(NMesh, 1:sps%emax, 0:sps%lmax))
    rnl(:,:,:) = 0.d0
    do n = 1, sps%emax
        do l = 0, min(n-1,sps%lmax)
        do i = 1, NMesh
          rnl(i,n,l) = hydrogen_radial_wf_norm(n, l, 1.d0, rmesh(i))
        end do
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
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
    deallocate(rnl)
    deallocate(rmesh, rwmesh)
  end subroutine set_coulomb_hydrogen

  subroutine set_kinetic_correction_laguerre(this, sps, zeta, NMesh, pmax)
    use AtLibrary, only: hc, m_e, gauss_legendre, &
        & Mom_laguerre_radial_wf_norm, alpha
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta, pmax
    integer, intent(in) :: NMesh
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    integer :: n, l, i
    real(8) :: r

    call gauss_legendre(0.d0, pmax, pmesh, pwmesh, NMesh)
    allocate(rnl_mom(NMesh, 0:sps%emax, 0:sps%lmax))
    rnl_mom(:,:,:) = 0.d0
    do n = 0, sps%emax
      do l = 0, sps%lmax
        do i = 1, NMesh
          rnl_mom(i,n,l) = Mom_laguerre_radial_wf_norm(n, l, 1.d0, pmesh(i))
        end do
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
          r = r + pwmesh(i) * rnl_mom(i,oa%n,oa%l) * rnl_mom(i,ob%n,ob%l) * &
              & (pmesh(i)*zeta)**4 * 0.125d0 / alpha**2
        end do
        this%m(a,b) = -r
        this%m(b,a) = -r
      end do
    end do

    deallocate(rnl_mom)
    deallocate(pmesh, pwmesh)
  end subroutine set_kinetic_correction_laguerre

  subroutine set_Darwin_laguerre(this, sps, zeta)
    use AtLibrary, only: laguerre_radial_wf, alpha
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    real(8) :: r

    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        if(oa%l /= 0) cycle
        r = 0.125d0 * laguerre_radial_wf(oa%n,0,1.d0/zeta,0.d0) * laguerre_radial_wf(ob%n,0,1.d0/zeta,0.d0) / alpha**2
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
  end subroutine set_Darwin_laguerre

  subroutine set_spin_orbit_laguerre(this, sps, zeta, NMesh, rmax)
    use AtLibrary, only: hc, m_e, gauss_legendre, &
        & fixed_point_quadrature, ln_gamma, laguerre, &
        & laguerre_radial_wf_norm, alpha, g_s
    type(OneBodyOperator), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    real(8), intent(in) :: zeta, rmax
    integer, intent(in) :: NMesh
    integer :: a, b
    type(EleSingleParticleOrbit), pointer :: oa, ob
    integer :: n, l, i
    real(8) :: ls, r

#ifdef gauss_laguerre
    call fixed_point_quadrature("laguerre", NMesh, rmesh, rwmesh, weight_renorm=.false., &
        & a_in=0.d0, b_in=2.d0, alpha_in=0.d0)
#else
    call gauss_legendre(0.d0, rmax, rmesh, rwmesh, NMesh)
#endif /* gauss_laguerre */
    allocate(rnl(NMesh, 0:sps%emax, 0:sps%lmax))
    rnl(:,:,:) = 0.d0
    do n = 0, sps%emax
      do l = 0, sps%lmax
        do i = 1, NMesh
#ifdef gauss_laguerre
          rnl(i,n,l) = exp( 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+2*l+3))) * &
              & laguerre(n,dble(2*l+2),2.d0*rmesh(i)) * (2.d0*rmesh(i))**(l+1) * sqrt(2.d0)
#else
          rnl(i,n,l) = laguerre_radial_wf_norm(n, l, 1.d0, rmesh(i))
#endif /* gauss_laguerre */
        end do
      end do
    end do

    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if(oa%j /= ob%j) cycle
        if(oa%l /= ob%l) cycle
        if(oa%l == 0) cycle
        ls = dble(oa%j)*0.5d0*(dble(oa%j)*0.5d0+1.d0) - dble(oa%l*(oa%l+1)) - 0.75d0
        r = 0.d0
        do i = 1, NMesh
          r = r + rnl(i,oa%n,oa%l) * rnl(i,ob%n,ob%l) * rwmesh(i) * zeta**3/ rmesh(i)**3
        end do
        r = r * ls * 0.125d0 * g_s / alpha**2
        this%m(a,b) = r
        this%m(b,a) = r
      end do
    end do
    deallocate(rnl)
    deallocate(rmesh, rwmesh)
  end subroutine set_spin_orbit_laguerre

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
end module OneBodyTerms
