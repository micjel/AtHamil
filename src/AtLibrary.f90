module AtLibrary
  use, intrinsic :: iso_c_binding
  implicit none

  real(8), parameter, public :: pi = 3.141592741012573d0 ! \pi
  real(8), parameter, public :: hc = 197.32705d0         ! \hbar c [eV nm]
  real(8), parameter, public :: m_e = 510.9989461 ! electron mass [keV]
  real(8), parameter, public :: alpha = 137.035999d0     ! electric fine structure constant
  real(8), parameter, public :: g_e =-2.002319304362d0  ! electron g-factor
  !
  ! C interfaces
  !
  interface
    ! 3-j symbol
    function coupling_3j(j1,j2,j3,m1,m2,m3) bind(c,name='gsl_sf_coupling_3j')
      import c_int, c_double
      real(c_double) :: coupling_3j
      integer(c_int), value, intent(in) :: j1,j2,j3,m1,m2,m3
    end function coupling_3j

    ! 6-j symbol
    function coupling_6j(j1,j2,j3,j4,j5,j6) bind(c,name='gsl_sf_coupling_6j')
      import c_int, c_double
      real(c_double) :: coupling_6j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6
    end function coupling_6j

    ! 9-j symbol
    function coupling_9j(j1,j2,j3,j4,j5,j6,j7,j8,j9) bind(c,name='gsl_sf_coupling_9j')
      import c_int, c_double
      real(c_double) :: coupling_9j
      integer(c_int), value, intent(in) :: j1,j2,j3,j4,j5,j6,j7,j8,j9
    end function coupling_9j

    ! factorial n!
    function factorial(n) bind(c,name='gsl_sf_fact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: factorial
    end function factorial

    ! double factorial n!!
    function double_factorial(n) bind(c,name='gsl_sf_doublefact')
      import c_double, c_int
      integer(c_int), value, intent(in) :: n
      real(c_double) :: double_factorial
    end function double_factorial

    ! Gamma function \Gamma(x)
    function gamma_function(x) bind(c,name='gsl_sf_gamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: gamma_function
    end function gamma_function

    ! log of Gamma function ln \Gamma(x)
    function ln_gamma(x) bind(c,name='gsl_sf_lngamma')
      import c_double
      real(c_double), value, intent(in) :: x
      real(c_double) :: ln_gamma
    end function ln_gamma

    ! spherical bessel function j_l(x)
    function spherical_bessel(l,x) bind(c,name='gsl_sf_bessel_jl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: spherical_bessel
    end function spherical_bessel

    ! Legendre polynomial P_l(x)
    function legendre_polynomial(l,x) bind(c,name='gsl_sf_legendre_Pl')
      import c_int, c_double
      integer(c_int), value, intent(in) :: l
      real(c_double), value, intent(in) :: x
      real(c_double) :: legendre_polynomial
    end function legendre_polynomial

    ! associated Laguerre polynomial L^{(a)}_{n}(x)
    function laguerre(n,a,x) bind(c,name='gsl_sf_laguerre_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, x
      real(c_double) :: laguerre
    end function laguerre

    ! Gegenbauer polynomial C^{(lambda)}_{n}(x)
    function Gegenbauer_polynomial(n,lambda,x) bind(c,name='gsl_sf_gegenpoly_n')
      import c_int, c_double
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: lambda, x
      real(c_double) :: gegenbauer_polynomial
    end function Gegenbauer_polynomial

    ! Gauss-Legendre quadrature
    function gauss_legendre_allocate(n) &
          & bind(c,name='gsl_integration_glfixed_table_alloc')
      import c_int, c_ptr
      integer(c_int), value, intent(in) :: n
      type(c_ptr) :: gauss_legendre_allocate
    end function gauss_legendre_allocate
    function gauss_legendre_ith_point_weight(a,b,i,xi,wi,t) &
          & bind(c,name='gsl_integration_glfixed_point')
      import c_int, c_double, c_ptr
      real(c_double), value, intent(in) :: a, b
      integer(c_int), value, intent(in) :: i
      real(c_double) :: xi, wi
      type(c_ptr), value, intent(in) :: t
      integer(c_int) :: gauss_legendre_ith_point_weight
    end function gauss_legendre_ith_point_weight
    subroutine gauss_legendre_release(t) &
          & bind(c,name='gsl_integration_glfixed_table_free')
      import c_ptr
      type(c_ptr), value :: t
    end subroutine gauss_legendre_release

    ! fixed-point quadrature
    function integration_fixed_allocate(T, n, a, b, alpha, beta) &
          & bind(c,name='gsl_integration_fixed_alloc')
      import c_int, c_double, c_ptr
      integer(c_int), value, intent(in) :: n
      real(c_double), value, intent(in) :: a, b, alpha, beta
      type(c_ptr), value :: T
      type(c_ptr) :: integration_fixed_allocate
    end function integration_fixed_allocate
    subroutine integration_fixed_free( workspace ) bind(c, name="gsl_integration_fixed_free")
      import c_ptr
      type(c_ptr), value :: workspace
    end subroutine integration_fixed_free
    function integration_fixed_n( workspace ) bind(c, name="gsl_integration_fixed_n")
      import c_ptr, c_int
      type(c_ptr), value :: workspace
      integer(c_int) :: integration_fixed_n
    end function integration_fixed_n
    function integration_fixed_nodes( workspace ) bind(c, name="gsl_integration_fixed_nodes")
      import c_ptr
      type(c_ptr), value :: workspace
      type(c_ptr) :: integration_fixed_nodes
    end function integration_fixed_nodes
    function integration_fixed_weights( workspace ) bind(c, name="gsl_integration_fixed_weights")
      import c_ptr
      type(c_ptr), value :: workspace
      type(c_ptr) :: integration_fixed_weights
    end function integration_fixed_weights

    ! open, read, write, and close gzip file (additional interface gzip_open below)
    ! When you pass string to C, you need add NULL (achar(0)) in the end of strings.
    ! It is done in gzip_open.
    function gz_open(filename, mode) bind(c, name='gzopen')
      import c_char, c_ptr
      character(c_char) :: filename(*), mode(*)
      type(c_ptr) :: gz_open
    end function gz_open
    function gzip_read(f, buf, len) bind(c, name='gzgets')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_read
    end function gzip_read
    function gzip_write(f, buf, len) bind(c, name='gzwrite')
      import c_int, c_char, c_ptr
      type(c_ptr), value :: f
      character(c_char) :: buf(*)
      integer(c_int), value, intent(in) :: len
      type(c_ptr) :: gzip_write
    end function gzip_write
    function gzip_close( f ) bind(c, name='gzclose')
      import c_ptr
      type(c_ptr), value :: f
      type(c_ptr) :: gzip_close
    end function gzip_close
  end interface
  type(c_ptr) :: integration_fixed_legendre; bind(c, name="gsl_integration_fixed_legendre") :: integration_fixed_legendre
  type(c_ptr) :: integration_fixed_chebyshev; bind(c, name="gsl_integration_fixed_chebyshev") :: integration_fixed_chebyshev
  type(c_ptr) :: integration_fixed_gegenbauer; bind(c, name="gsl_integration_fixed_gegenbauer") :: integration_fixed_gegenbauer
  type(c_ptr) :: integration_fixed_jacobi; bind(c, name="gsl_integration_fixed_jacobi") :: integration_fixed_jacobi
  type(c_ptr) :: integration_fixed_laguerre; bind(c, name="gsl_integration_fixed_laguerre") :: integration_fixed_laguerre
  type(c_ptr) :: integration_fixed_hermite; bind(c, name="gsl_integration_fixed_hermite") :: integration_fixed_hermite
  type(c_ptr) :: integration_fixed_exponential; bind(c, name="gsl_integration_fixed_exponential") :: integration_fixed_exponential
  type(c_ptr) :: integration_fixed_rational; bind(c, name="gsl_integration_fixed_rational") :: integration_fixed_rational
  type(c_ptr) :: integration_fixed_chebyshev2; bind(c, name="gsl_integration_fixed_chebyshev2") :: integration_fixed_chebyshev2
  !
  ! end C interfaces
  !
  interface str
    procedure :: i2str
    procedure :: real2str
    procedure :: str2str
  end interface str
  integer, private, parameter :: Lenc = 256
contains

  subroutine skip_comment(nfile, comment)
    integer,intent(in)::nfile
    character(*), intent(in) :: comment
    character(20) :: line
    read(nfile,'(a)') line
    do while (find(line, comment))
      read(nfile,'(a)') line
    end do
    backspace(nfile)
  end subroutine skip_comment

  pure real(8) function delta(i1, i2)
    integer, intent(in) :: i1, i2
    delta = 0.d0
    if(i1 == i2) delta = 1.d0
  end function delta

  pure function hat(j) result(r)
    real(8) :: r
    integer, intent(in) :: j
    r = dsqrt(dble(j + 1))
  end function hat

  ! 3-j symbol for Wigner-Eckert theorem
  real(8) function geometry_part(jbra, jop, jket, mbra, mop, mket) result(r)
    integer, intent(in) :: jbra, mbra, jop, mop, jket, mket
    r = (-1.d0) ** ((jbra - mbra)/2) * &
      & tjs(jbra, jop, jket, -mbra, mop, mket)
  end function geometry_part

  pure real(8) function red_r_l(n1, l1, n2, l2) result(rl)
    integer, intent(in) :: n1, l1, n2, l2
    if (n1 == n2 .and. l1 == l2-1) then
      rl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
    elseif (n1 == n2-1 .and. l1 == l2+1) then
      rl = -dsqrt(dble(l2 + 1)*dble(n2))
    elseif (n1 == n2+1 .and. l1 == l2-1) then
      rl = dsqrt(dble(l2)*(dble(n2  + 1)))
    elseif (n1 == n2 .and. l1==l2+1) then
      rl = dsqrt(dble(l2+1)*(dble(n2 +l2)+1.5d0))
    else
      rl = 0.d0
    end if
  end function red_r_l

  pure real(8) function red_nab_l(n1, l1, n2, l2) result(nl)
    integer, intent(in) :: n1, l1, n2, l2
    if(n1 == n2 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*(dble(n2 + l2) + 1.5d0))
    elseif(n1 == n2-1 .and. l1 == l2+1) then
      nl = -dsqrt(dble(l2 + 1)*dble(n2))
    elseif(n1 == n2 .and. l1 == l2-1) then
      nl = -dsqrt(dble(l2)*(dble(n2 + l2) + 0.5d0))
    elseif(n1 == n2+1 .and. l1==l2-1) then
      nl = -dsqrt(dble(l2)*dble(n2 + 1))
    else
      nl = 0.d0
    end if
  end function red_nab_l

  pure logical function triag(i,j,k)
      implicit none
      integer,intent(in)::i,j,k
      triag = ((i-(j+k))*(i-abs(j-k)) > 0)
  end function triag

  function dcg(j1, m1, j2, m2, j3, m3) result(s)
    !
    !  Clebsch-Gordan coefficient
    !  dcg(j1, m1, j2, m2, j3, m3)
    !  = ((j1)/2, (m1)/2, (j2)/2, (m2)/2 | (j3)/2, (m3)/2)
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    real(8) :: s
    s = coupling_3j(j1,j2,j3,m1,m2,-m3) * hat(j3) * (-1.d0) ** ((j1-j2+m3)/2)
  end function dcg

  function tjs(j1, j2, j3, m1, m2, m3) result(r)
    real(8) :: r
    integer, intent(in) :: j1, j2, j3, m1, m2, m3
    r = coupling_3j(j1,j2,j3,m1,m2,m3)
  end function tjs

  function sjs(j1, j2, j3, l1, l2, l3) result(s)
    !
    !  6j coefficient
    !  d6j(j1, j2, j3, l1, l2, l3) = {(j1)/2 (j2)/2 (j3)/2}
    !                                {(l1)/2 (l2)/3 (l3)/2}
    integer, intent(in) :: j1, j2, j3, l1, l2, l3
    real(8) :: s
    s = coupling_6j(j1,j2,j3,l1,l2,l3)
  end function sjs

  function snj(j11, j12, j13, j21, j22, j23, j31, j32, j33) result(s)
    !
    !  9j coefficient
    !  d9j(j11, j12, j13, j21, j22, j23, j31, j32, j33)
    !    {(j11)/2 (j12)/2 (j13)/2}
    !  = {(j21)/2 (j22)/2 (j23)/2}
    !    {(j31)/2 (j32)/2 (j33)/2}
    integer, intent(in) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
    real(8) :: s
    s = coupling_9j(j11,j12,j13,j21,j22,j23,j31,j32,j33)
  end function snj

  function ho_radial_wf(n,l,anu,r) result(s)
    ! R = sqrt( 2 * nu**3 * Gamma(n+1) / Gamma(n+l+1.5) ) x^{l} e^{ -x^2/2 } L^{n}_{l+0.5}( x^2 )
    ! x = r * nu or p * nu
    ! nu = sqrt(m w / h) or sqrt( h / mw )
    ! Note that anu = nu^2
    integer,intent(in) :: n,l
    real(8),intent(in) :: anu,r
    real(8) :: s
    real(8) :: prefact, exp_component, nu, x, a
    nu = sqrt(anu)
    x = nu * r
    a = dble(l)+0.5d0
    prefact = sqrt( 2.d0 * nu**3 )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+l)+1.5d0) - 0.5d0*x**2 &
        & + dble(l) * log(x)
    s = prefact * exp(exp_component) * laguerre(n,a,x**2)
  end function ho_radial_wf

  function ho_radial_wf_norm(n,l,anu,r) result(s)
    ! R = sqrt( 2 * nu * Gamma(n+1) / Gamma(n+l+1.5) ) x^{l+1} e^{ -x^2/2 } L^{n}_{l+0.5}( x^2 )
    ! x = r * nu or p * nu
    ! nu = sqrt(m w / h) or sqrt( h / mw )
    ! Note that anu = nu^2
    integer,intent(in) :: n,l
    real(8),intent(in) :: anu,r
    real(8) :: s
    real(8) :: prefact, exp_component, nu, x, a
    nu = sqrt(anu)
    x = nu * r
    a = dble(l)+0.5d0
    prefact = sqrt( 2.d0 * nu )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - 0.5d0*ln_gamma(dble(n+l)+1.5d0) - 0.5d0*x**2 &
        & + dble(l+1) * log(x)
    s = prefact * exp(exp_component) * laguerre(n,a,x**2)
  end function ho_radial_wf_norm

  function hydrogen_radial_wf(n,l,a_star,r) result(s)
    ! R_{nl} = sqrt( (2/n a* )**3 (n-l-1)!/2n(n+l)! ) x^{l} e^{ -x/2 } L^{2l+1}_{n-l-1}( x )
    ! x = 2 * r / (n*a*)
    ! a* is reduced Bohr radius
    integer,intent(in) :: n,l
    real(8),intent(in) :: a_star,r
    real(8) :: s
    real(8) :: prefact, exp_component,  x
    x = 2.d0 * r / ( a_star * dble(n))
    prefact = sqrt( 4.d0/ a_star**3 ) / dble(n**2)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) - 0.5d0*x + dble(l)*log(x)
    s = prefact * exp(exp_component) * laguerre(n-l-1,dble(2*l+1),x)
  end function hydrogen_radial_wf

  function hydrogen_radial_wf_norm(n,l,a_star,r) result(s)
    ! R_{nl} = sqrt( (2/n a* ) (n-l-1)!/2n(n+l)! ) x^{l+1} e^{ -x/2 } L^{2l+1}_{n-l-1}( x )
    ! x = 2 * r / (n*a*)
    ! a* is reduced Bohr radius
    integer,intent(in) :: n,l
    real(8),intent(in) :: a_star,r
    real(8) :: s
    real(8) :: prefact, exp_component,  x
    x = 2.d0 * r / ( a_star * dble(n))
    prefact = sqrt( 1.d0/a_star ) / dble(n)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) - 0.5d0*x + dble(l+1)*log(x)
    s = prefact * exp(exp_component) * laguerre(n-l-1,dble(2*l+1),x)
  end function hydrogen_radial_wf_norm

  function Laguerre_radial_wf(n,l,b,r) result(s)
    ! From Eq.(13) in A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).
    ! b is length scale [L]
    integer,intent(in) :: n,l
    real(8),intent(in) :: b,r
    real(8) :: s
    real(8) :: prefact, exp_component, x, x2
    x = r / b
    x2 = 2.d0 * x
    prefact = sqrt( 2.d0 / b ) * 2.d0 / b
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - &
        & 0.5d0*ln_gamma(dble(n+2*l+3)) - x
    if(abs(x2) > 1.d-16) then
      exp_component = exp_component  + dble(l)*log(x2)
    end if
    s = prefact * exp(exp_component) * laguerre(n,dble(2*l+2),x2)
  end function Laguerre_radial_wf

  function Laguerre_radial_wf_norm(n,l,b,r) result(s)
    ! From Eq.(13) in A. E. McCoy and M. A. Caprio, J. Math. Phys. 57, (2016).
    ! b is length scale [L]
    integer,intent(in) :: n,l
    real(8),intent(in) :: b,r
    real(8) :: s
    real(8) :: prefact, exp_component, x, x2
    x = r / b
    x2 = 2.d0 * x
    prefact = sqrt( 2.d0 / b )
    exp_component = 0.5d0*ln_gamma(dble(n+1)) - &
        & 0.5d0*ln_gamma(dble(n+2*l+3)) - x  + dble(l+1) * log(x2)
    s = prefact * exp(exp_component) * laguerre(n,dble(2*l+2),x2)
  end function Laguerre_radial_wf_norm

  function Mom_Laguerre_radial_wf_norm(n,l,b,p) result(s)
    ! b is length scale [L]
    integer,intent(in) :: n,l
    real(8),intent(in) :: b,p
    real(8) :: s
    real(8) :: exp_prefact, exp_component, x, xi
    integer :: k
    xi = b * p * 0.5d0
    x = 0.5d0 / sqrt(0.25d0 + xi*xi)
    s = 0.d0
    exp_prefact = 0.5d0*ln_gamma(dble(n+1)) + 0.5d0*ln_gamma(dble(n+2*l+3)) + &
      & dble(2*l+3)*log(2.d0*x) + dble(l)*log(2.d0*xi) + ln_gamma(dble(l+1))
    do k = 0, n
      exp_component = dble(k)*log(2.d0*x) - ln_gamma(dble(n+1-k)) - ln_gamma(dble(2*l+k+3)) + exp_prefact
      s = s + (-1.d0)**k * dble(k+1) * Gegenbauer_polynomial(k+1,dble(l+1),x) * exp(exp_component)
    end do
    s = s * sqrt(b / pi) * xi
  end function Mom_Laguerre_radial_wf_norm

  function hydrogen_radial_wf_mom_norm(n,l,a_star,p) result(s)
    ! from wikipedia https://en.wikipedia.org/wiki/Hydrogen_atom
    ! Do not forget phase (-i)^{l} when you calclate the matrix elements!!
    ! x = p * (n*a*)
    ! a* is reduced Bohr radius
    integer,intent(in) :: n, l
    real(8),intent(in) :: a_star, p
    real(8) :: s
    real(8) :: prefact, exp_component, x, z
    x = p * dble(n) * a_star
    z = (x**2-1.d0) / (x**2+1.d0)
    prefact = sqrt( 2.d0/pi * a_star) * dble(n)
    exp_component = 0.5d0*ln_gamma(dble(n-l)) - 0.5d0*ln_gamma(dble(n+l+1)) + dble(2*l+2)*log(2.d0) + &
        & ln_gamma(dble(l+1)) + dble(l+1)*log(x) - dble(l+2)*log(x**2+1)
    s = prefact * exp(exp_component) * Gegenbauer_polynomial(n-l-1,dble(l+1),z)
  end function hydrogen_radial_wf_mom_norm

  subroutine gauss_legendre(x1,x2,x,w,n)
    integer, intent(in) :: n
    real(8), intent(in) :: x1, x2
    real(8), intent(out), allocatable :: x(:), w(:)
    real(8) :: xi, wi
    integer :: info, i
    type(c_ptr) :: t

    if(allocated(x)) deallocate(x)
    if(allocated(w)) deallocate(w)
    allocate(x(n))
    allocate(w(n))
    t = gauss_legendre_allocate(n)
    do i = 1, n
      info = gauss_legendre_ith_point_weight(x1,x2,i-1,xi,wi,t)
      x(i) = xi
      w(i) = wi
    end do
    call gauss_legendre_release(t)
  end subroutine gauss_legendre

  subroutine fixed_point_quadrature(quad_name, n, x, w, a_in, b_in, alpha_in, beta_in, weight_renorm)
    character(*), intent(in) :: quad_name
    integer, intent(in) :: n
    real(8), intent(in), optional :: a_in, b_in, alpha_in, beta_in
    logical, intent(in), optional :: weight_renorm
    real(8), allocatable :: x(:), w(:)
    real(8) :: a=0.d0, b=0.d0, alpha=0.d0, beta=0.d0
    integer :: i
    type(c_ptr) :: workspace, nodes, weights
    real(c_double), pointer :: x_(:), w_(:)

    if(allocated(x)) deallocate(x)
    if(allocated(w)) deallocate(w)
    allocate(x(n))
    allocate(w(n))
    if(present(a_in)) a = a_in
    if(present(b_in)) b = b_in
    if(present(alpha_in)) alpha = alpha_in
    if(present(beta_in)) beta = beta_in

    select case(quad_name)
    case("legendre","Legendre")
      workspace = integration_fixed_allocate(integration_fixed_legendre, n, a, b, alpha, beta)
    case("chebyshev1","Chebyshev1","Chebyshev type 1", "Chebyshev Type 1")
      workspace = integration_fixed_allocate(integration_fixed_chebyshev, n, a, b, alpha, beta)
    case("gegenbauer","Gegenbauer")
      workspace = integration_fixed_allocate(integration_fixed_gegenbauer, n, a, b, alpha, beta)
    case("jacobi","Jacobi")
      workspace = integration_fixed_allocate(integration_fixed_jacobi, n, a, b, alpha, beta)
    case("laguerre","Laguerre")
      workspace = integration_fixed_allocate(integration_fixed_laguerre, n, a, b, alpha, beta)
    case("hermite","Hermite")
      workspace = integration_fixed_allocate(integration_fixed_hermite, n, a, b, alpha, beta)
    case("exponential","Exponential")
      workspace = integration_fixed_allocate(integration_fixed_exponential, n, a, b, alpha, beta)
    case("rational","Rational")
      workspace = integration_fixed_allocate(integration_fixed_rational, n, a, b, alpha, beta)
    case("chebyshev2","Chebyshev2","Chebyshev type 2", "Chebyshev Type 2")
      workspace = integration_fixed_allocate(integration_fixed_chebyshev2, n, a, b, alpha, beta)
    case default
      write(*,"(a)") "Unknown quadrature name"
      stop
    end select

    nodes = integration_fixed_nodes( workspace )
    weights = integration_fixed_weights( workspace )
    call c_f_pointer(nodes, x_, [n])
    call c_f_pointer(weights, w_, [n])
    x(:) = x_(:)
    w(:) = w_(:)
    call integration_fixed_free(workspace)
    if(.not. present(weight_renorm)) return
    if(.not. weight_renorm) return
    do i = 1, n
      select case(quad_name)
      case("legendre","Legendre")
        w(i) = w(i) * 1.d0
      case("chebyshev1","Chebyshev1","Chebyshev type 1", "Chebyshev Type 1")
        w(i) = w(i) * sqrt((b-x(i)) * (x(i)-a))
      case("gegenbauer","Gegenbauer")
        w(i) = w(i) / ( (b-x(i)) * (x(i)-a) )**alpha
      case("jacobi","Jacobi")
        w(i) = w(i) / ( (b-x(i))**alpha * (x(i)-a)**beta )
      case("laguerre","Laguerre")
        w(i) = w(i) * exp( b * (x(i)-a) ) / (x(i)-a)**alpha
      case("hermite","Hermite")
        w(i) = w(i) * exp( b * (x(i)-a)**2 ) / abs(x(i)-a)**alpha
      case("exponential","Exponential")
        w(i) = w(i) / abs( x(i)-(a+b)*0.5d0 )**alpha
      case("rational","Rational")
        w(i) = w(i) / ( (x(i)-a)**alpha * (x(i)+b)*beta )
      case("chebyshev2","Chebyshev2","Chebyshev type 2", "Chebyshev Type 2")
        w(i) = w(i) / sqrt( (b-x(i)) * (x(i)-a) )
      case default
        write(*,"(a)") "Unknown quadrature name"
        stop
      end select
    end do
  end subroutine fixed_point_quadrature

  function gzip_open( filename, mode ) result(p)
    character(*), intent(in) :: filename, mode
    type(c_ptr) :: p
    p = gz_open(trim(filename)//achar(0), trim(mode)//achar(0))
  end function gzip_open

  function gzip_writeline( f, buf, len) result(p)
    type(c_ptr) :: f, p
    character(*), intent(in) :: buf
    integer, intent(in) :: len
    p = gzip_write( f, trim(buf)//achar(10), len+1)
  end function gzip_writeline

  function gzip_readline( f, buf, len) result(p)
    ! note
    ! len_trim returns length removed space (32 in ascii code)
    ! -2 means removing null (0) and line feed (10)
    type(c_ptr) :: f, p
    character(*), intent(inout) :: buf
    integer, intent(in) :: len
    p = gzip_read( f, buf, len)
    buf = buf(1:len_trim(buf) - 2)
  end function gzip_readline

  logical function find(str, key) result(r)
    character(*), intent(in) :: str, key
    integer :: i
    i = index(str, key)
    if(i == 0) r = .false.
    if(i /= 0) r = .true.
  end function find

  subroutine mkdir(dir)
    character(*), intent(in) :: dir
    character(Lenc) :: comm
    integer :: is, ie, idxdir
    integer :: lnblnk

    idxdir = lnblnk(dir)
    is = 1; ie = is + len('if [ ! -d ')
    write(comm(is:ie), '(a)') 'if [ ! -d '
    is = ie + 1; ie = is + idxdir
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1
    ie = is + len(' ]; then (echo "Creating directory, ')
    write(comm(is:ie), '(a)') ' ]; then (echo "Creating directory, '
    is = ie + 1; ie = is + idxdir
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1; ie = is + len('"; mkdir -p ') - 1
    write(comm(is:ie), '(a)') '"; mkdir -p '
    is = ie + 1; ie = is + idxdir - 1
    write(comm(is:ie), '(a)') dir(1:idxdir)
    is = ie + 1; ie = is + 4
    write(comm(is:ie), '(a)') '); fi'
    !write(*,'(a)') comm(1:ie)
    call system(comm(1:ie))
  end subroutine mkdir

  function i2str(i)
    character(:), allocatable :: i2str
    character(Lenc) :: ist
    integer, intent(in) :: i
    write(ist,*) i
    i2str=adjustl(trim(ist))
  end function i2str

  function real2str(r) result(str)
    character(:), allocatable :: str
    character(Lenc) ::  rst
    real(8), intent(in) :: r
    integer :: l
    l = int(log10(abs(r) + 1.d-8))
    if(l >= 1) then
      write(rst, *) int(r)
    elseif(l < 1 .and. l >= 0) then
      write(rst, '(f10.2)') r
    elseif(l < 0 .and. l >= -1) then
      write(rst, '(f10.3)') r
    elseif(l < -1 .and. l >= -2) then
      write(rst, '(f10.4)') r
    elseif(l < -2 .and. l >= -3) then
      write(rst, '(f10.5)') r
    elseif(l < -3 .and. l >= -4) then
      write(rst, '(f10.6)') r
    elseif(l < -4 .and. l >= -5) then
      write(rst, '(f10.7)') r
    elseif(l < -5 .and. l >= -6) then
      write(rst, '(f10.8)') r
    elseif(l < -6 .and. l >= -7) then
      write(rst, '(f10.9)') r
    else
      rst = "0"
    end if
    str = adjustl(trim(rst))
  end function real2str

  function str2str(str)
    character(:), allocatable :: str2str
    character(*), intent(in) :: str
    str2str=adjustl(trim(str))
  end function str2str

  function isfile(f) result(ex)
    character(*), intent(in) :: f
    logical :: ex
    inquire(file=f, exist=ex)
  end function isfile
end module AtLibrary

