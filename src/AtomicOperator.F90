module AtomicOperator
  use omp_lib
  use LinAlgLib
  use ElectronTwoBodySpace
  use EleSingleParticleState
  use OneBodyTerms
  use TwoBodyTerms
  implicit none
  type :: Op
    type(OneBodyOperator) :: one
    type(TwoBodyOperator) :: two
    type(EleOrbits), pointer :: sps
    type(EleTwoBodySpace), pointer :: ms
    integer :: jrank = 0, prank = 1
    character(:), allocatable :: OpName
  contains
    procedure :: InitOp
    procedure :: FinOp
    procedure :: SetOp
    procedure :: WriteOp
    generic :: init => InitOp
    generic :: fin => FinOp
    generic :: set => SetOp
  end type Op

  real(8), private, allocatable :: sp_radial_integrals(:,:,:,:)
contains

  subroutine InitOp( this, OpName, two, jrank, prank )
    class(Op), intent(inout) :: this
    character(*), intent(in) :: OpName
    type(EleTwoBodySpace), intent(in), target :: two
    integer, intent(in), optional :: jrank, prank
    this%ms => two
    this%sps => two%sps
    this%OpName = OpName
    this%jrank = jrank
    this%prank = prank
    call this%one%init(two%sps, jrank, prank)
    call this%two%init(two, jrank, prank)
  end subroutine InitOp

  subroutine FinOp( this )
    class(Op), intent(inout) :: this
    call this%one%fin()
    call this%two%fin()
    this%sps => null()
    this%ms => null()
    this%jrank = 0
    this%prank = 0
    deallocate(this%OpName)
  end subroutine FinOp

  subroutine SetOp( this, NMesh, pmax )
    class(Op), intent(inout) :: this
    integer, intent(in), optional :: NMesh
    real(8), intent(in), optional :: pmax
    select case( this%OpName )
    case("MassShift")
      call set_mass_shift( this, NMesh, pmax )
    case("FieldShift")
      call set_field_shift( this )
    end select
  end subroutine SetOp

  subroutine WriteOp( this, f )
    use AtLibrary, only: find
    class(Op), intent(inout) :: this
    character(*), intent(in) :: f
    if(find(f,".snt")) then
      call write_operator_snt( this, f )
      return
    end if

    if(find(f,".op.gz")) then
      call write_operator_gzip( this, f )
      return
    end if
  end subroutine WriteOp

  subroutine set_mass_shift( this, NMesh, pmax )
    class(Op), intent(inout) :: this
    integer, intent(in) :: NMesh
    real(8), intent(in) :: pmax

    call set_psquare_one_body( this%one, this%ms%zeta, NMesh, pmax )
    call set_p_dot_p_two_body( this%two, NMesh, pmax )
  end subroutine set_mass_shift

  subroutine set_field_shift( this )
    class(Op), intent(inout) :: this
    call set_field_shift_0th_order( this%one, this%ms%zeta )
  end subroutine set_field_shift

  subroutine set_psquare_one_body( this, zeta, NMesh, pmax )
    use AtLibrary, only: Mom_laguerre_radial_wf_norm
    type(OneBodyOperator), intent(inout) :: this
    real(8), intent(in) :: zeta, pmax
    integer, intent(in) :: NMesh
    integer :: na, nb, l, i
    integer :: a, b
    type(EleOrbits), pointer :: sps
    type(EleSingleParticleOrbit), pointer :: oa, ob
    real(8) :: r
    real(8), allocatable :: Mesh(:), WMesh(:)
    sps => this%sps
    allocate( sp_radial_integrals( 0:sps%emax, 0:sps%emax, 0:sps%emax, 0:sps%emax))
    sp_radial_integrals(:,:,:,:) = 0.d0
    call gauss_legendre(0.d0, pmax, Mesh, WMesh, NMesh)
    do na = 0, sps%emax
      do nb = 0, sps%emax
        do l = 0, sps%emax
          r = 0.d0
          do i = 1, NMesh
            r = r + 0.5d0 * (Mesh(i)*zeta)**2 * WMesh(i) * &
                & Mom_laguerre_radial_wf_norm(na,l,1.d0,Mesh(i)) * &
                & Mom_laguerre_radial_wf_norm(nb,l,1.d0,Mesh(i))
          end do
          sp_radial_integrals(na,nb,l,l) = r
        end do
      end do
    end do
    deallocate( Mesh, WMesh )

    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if( oa%l /= ob%l ) cycle
        if( oa%j /= ob%j ) cycle
        this%m(a,b) = sp_radial_integrals(oa%n,ob%n,oa%l,oa%l)
        this%m(b,a) = sp_radial_integrals(ob%n,oa%n,oa%l,oa%l)
      end do
    end do
    deallocate(sp_radial_integrals)
  end subroutine set_psquare_one_body

  subroutine set_p_dot_p_two_body( this, NMesh, pmax )
    use AtLibrary, only: laguerre_radial_wf_norm, d_laguerre_radial_wf_norm, &
        & Mom_laguerre_radial_wf_norm
    type(TwoBodyOperator), intent(inout) :: this
    real(8), intent(in) :: pmax
    integer, intent(in) :: NMesh
    integer :: na, nb, la, lb, i, ch, J
    integer :: bra, ket, a, b, c, d
    type(EleOrbits), pointer :: sps
    type(EleTwoBodySpace), pointer :: ms
    type(EleTwoBodyChan), pointer :: tbc
    real(8) :: rr, r, fact, zeta
    real(8), allocatable :: Mesh(:), WMesh(:)
    ms => this%ms
    sps => this%sps
    zeta = ms%zeta
    allocate( sp_radial_integrals( 0:sps%emax, 0:sps%emax, 0:sps%emax, 0:sps%emax) )
    sp_radial_integrals(:,:,:,:) = 0.d0
    call gauss_legendre(0.d0, pmax, Mesh, WMesh, NMesh)
    do na = 0, sps%emax
      do nb = 0, sps%emax
        do la = 0, sps%emax
          do lb = 0, sps%emax
            fact = 0.d0
            if( la == lb+1 ) fact = sqrt(dble(lb+1))
            if( la == lb-1 ) fact =-sqrt(dble(lb))
            r = 0.d0
            rr = 0.d0
            if( la == lb+1 ) then
              do i = 1, NMesh
                rr = rr + zeta * WMesh(i) * Laguerre_radial_wf_norm(na,la,1.d0,Mesh(i)) * &
                    & (                 d_Laguerre_radial_wf_norm(nb,lb,1.d0,Mesh(i) ) - &
                    & dble(lb) / Mesh(i) * Laguerre_radial_wf_norm( nb,lb,1.d0,Mesh(i)) ) * fact
              end do
            end if
            if( la == lb-1 ) then
              do i = 1, NMesh
                rr = rr + zeta * WMesh(i) * Laguerre_radial_wf_norm(na,la,1.d0,Mesh(i)) * &
                    & (                    d_Laguerre_radial_wf_norm( nb,lb,1.d0,Mesh(i) ) + &
                    & dble(lb+1) / Mesh(i) * Laguerre_radial_wf_norm( nb,lb,1.d0,Mesh(i))) * fact
              end do
            end if
            !do i = 1, NMesh
            !  r = r + (Mesh(i)*zeta) * WMesh(i) * &
            !      & Mom_laguerre_radial_wf_norm(na,la,1.d0,Mesh(i)) * &
            !      & Mom_laguerre_radial_wf_norm(nb,lb,1.d0,Mesh(i)) * fact * (-1.d0)**( (la-lb+1)/2 )
            !end do
            !write(*,*) na,la,nb,lb,r, rr
            sp_radial_integrals(na,nb,la,lb) = rr
          end do
        end do
      end do
    end do
    deallocate( Mesh, WMesh )

    do ch = 1, ms%NChan
      tbc => ms%jp(ch)
      J = tbc%J
      do bra = 1, tbc%n_state
        a = tbc%n2label1(bra)
        b = tbc%n2label2(bra)
        do ket = 1, bra
          c = tbc%n2label1(ket)
          d = tbc%n2label2(ket)
          r = p_dot_p( a, b, c, d, J )
          this%MatCh(ch,ch)%m(bra,ket) = r
        end do
      end do
    end do
    deallocate(sp_radial_integrals)
  contains
    function p_dot_p( a, b, c, d, J ) result(r)
      integer, intent(in) :: a, b, c, d, J
      type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
      real(8) :: r, norm
      oa => sps%orb(a)
      ob => sps%orb(b)
      oc => sps%orb(c)
      od => sps%orb(d)
      norm = 1.d0
      if(a==b) norm = norm * dsqrt(0.5d0)
      if(c==d) norm = norm * dsqrt(0.5d0)
      r = - (-1.d0)**( (ob%j+oc%j)/2+J ) * sjs(oa%j, ob%j, 2*J, od%j, oc%j, 2) * red_nab_j(a,c) * red_nab_j(b,d) - &
          & (-1.d0)**( (ob%j+oc%j)/2 ) * sjs(oa%j, ob%j, 2*J, oc%j, od%j, 2) * red_nab_j(a,d) * red_nab_j(b,c)
    end function p_dot_p

    function red_nab_j(a, b) result(r)
      integer, intent(in) :: a, b
      real(8) :: r
      type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
      oa => sps%orb(a)
      ob => sps%orb(b)
      oc => sps%orb(c)
      od => sps%orb(d)
      r = (-1.d0)**((3+2*oa%l+ob%j)/2) * &
          & dsqrt(dble(oa%j + 1) * dble(ob%j + 1)) * &
          & sjs(oa%j, 2, ob%j, 2*ob%l, 1, 2*oa%l) * &
          & sp_radial_integrals(oa%n,oa%l,ob%n,ob%l)
    end function red_nab_j
  end subroutine set_p_dot_p_two_body

  subroutine set_field_shift_0th_order( this, zeta )
    type(OneBodyOperator), intent(inout) :: this
    real(8), intent(in) :: zeta
    integer :: a
    type(EleSingleParticleOrbit), pointer :: oa
    real(8) :: r
    type(EleOrbits), pointer :: sps
    sps => this%sps
    do a = 1, sps%norbs
      oa => sps%orb(a)
      if(oa%l /= 0 ) cycle
      r = 2.d0 * pi * &
          & laguerre_radial_wf(oa%n,oa%l,1.d0/zeta,0.d0) * &
          & laguerre_radial_wf(oa%n,oa%l,1.d0/zeta,0.d0) / 3.d0
      this%m(a,a) = r
    end do
  end subroutine set_field_shift_0th_order

  subroutine write_operator_snt( this, f )
    type(Op), intent(in) :: this
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
    write(wunit, "(3a)") "# Atomic operator with ", trim(sps%basis), " basis"
    write(wunit, "(a,f6.4,a)") "# length scale is a_0 / ", ms%zeta, " with bohr radius a_0"
    write(wunit, "(a,f8.4,a)") "# Corresponding hw is ", &
        & this%ms%zeta**2 * 27.211396, " eV"
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
        if(abs(this%one%m(bra,ket)) < 1.d-12) cycle
        cnt = cnt + 1
      end do
    end do

    write(wunit, '(a)') '####### one-body term'
    write(wunit, '(i5, i3)') cnt, 0
    write(wunit, '(a)') '### a b <a||b>'
    do bra = 1, sps%norbs
      do ket = 1, bra
        if(sps%orb(bra)%l /= sps%orb(ket)%l) cycle
        if(sps%orb(bra)%j /= sps%orb(ket)%j) cycle
        if(abs(this%one%m(bra,ket)) < 1.d-12) cycle
        write(wunit,'(2i5,1es20.10)') bra,ket, this%one%m(bra,ket)
      end do
    end do

    cnt = 0
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
          if( abs(this%two%GetTBME(a,b,c,d,J)) < 1.d-12) cycle
          cnt = cnt + 1
        end do
      end do
    end do

    write(wunit, '(a)') '##### two-body term #####'
    write(wunit, '(i15, i3)') cnt, 0
    write(wunit, '(a)') '## a, b, c, d, J, <ab||cd>'
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
          if( abs(this%two%GetTBME(a,b,c,d,J)) < 1.d-12) cycle
          write(wunit, '(5i5, 6es20.10)') a, b, c, d, J, &
              & this%two%GetTBME(a,b,c,d,J)
        end do
      end do
    end do
    close(wunit)
  end subroutine write_operator_snt

  subroutine write_operator_gzip( this, f )
    use AtLibrary, only: gzip_open, gzip_writeline, gzip_close
    type(Op), intent(in) :: this
    character(*), intent(in) :: f
    type(EleTwoBodySpace), pointer :: ms
    type(EleOrbits), pointer :: sps
    integer :: a, b, c, d, bra, ket, Jab, Jcd, ch
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: cnt, n, line
    character(256) :: header
    type(c_ptr) :: fp, err
    character(kind=c_char, len=200) :: buffer
    integer :: nbuf = 10000000, nrest
    character(12) :: cfmt
    real(8), allocatable :: v_buffer(:)
    ms => this%ms
    sps => ms%sps
#ifdef VERSION
    header = "# Generated by NuHamil (Tokyo code), "// trim(VERSION)
#else
    header = "# Generated by NuHamil (Tokyo code)"
#endif
    fp = gzip_open(f,"wt")
    err = gzip_writeline(fp, trim(header), len_trim(header))
    write(buffer,"(5i4)") this%jrank, this%prank, ms%emax, ms%e2max, sps%lmax
    err = gzip_writeline(fp, trim(buffer), len_trim(buffer))

    do bra = 1, sps%norbs
      do ket = 1, bra
        write(buffer,"(1es20.10)") this%one%m(bra,ket)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end do
    end do

    allocate( v_buffer(10*nbuf) )
    cnt = 0
    do a = 1, sps%norbs
      oa => sps%orb(a)
      do b = 1, a
        ob => sps%orb(b)
        if( oa%e + ob%e > ms%e2max ) cycle
        do c = 1, a
          oc => sps%orb(c)
          do d = 1, c
            od => sps%orb(d)
            if( oc%e + od%e > ms%e2max ) cycle
            if( (-1)**(oa%l+ob%l+oc%l+od%l)*this%prank==-1 ) cycle

            do Jab = abs(oa%j-ob%j)/2, (oa%j+ob%j)/2
              do Jcd = abs(oc%j-od%j)/2, (oc%j+od%j)/2
                if( triag( Jab, Jcd, this%jrank )) cycle
                cnt = cnt + 1
                v_buffer(cnt) = this%two%GetTBME(a,b,c,d,Jab,Jcd)
                if( cnt == 10*nbuf ) then
                  do line = 1, nbuf
                    write(buffer, "(10es20.10)") v_buffer(10*(line-1)+1:10*line)
                    err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
                  end do
                  cnt = 0
                end if
              end do
            end do
          end do
        end do
      end do
    end do

    if( cnt /= 0 ) then
      do line = 1, cnt/10
        write(buffer, "(10es20.10)") v_buffer(10*(line-1)+1:10*line)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end do

      nrest = cnt - (cnt/10) * 10
      if(nrest > 0) then
        cfmt = '(xes20.10)'
        write(cfmt(2:2),'(i1)') nrest
        write(buffer,cfmt) v_buffer((cnt/10)*10+1:cnt)
        err = gzip_writeline(fp, trim(buffer), len_trim(buffer))
      end if

    end if
    err = gzip_close(fp)
    deallocate(v_buffer)
  end subroutine write_operator_gzip
end module AtomicOperator
