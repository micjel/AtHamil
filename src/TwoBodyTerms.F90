module TwoBodyTerms
  use omp_lib
  use LinAlgLib
  use ElectronTwoBodySpace
  use EleSingleParticleState
  use AtLibrary
  implicit none

  type, extends(DMat) :: TwoBodyOpChan
    type(EleOrbits), pointer :: sps
    type(EleTwoBodyChan), pointer :: chbra, chket
    logical :: zero = .true.
  contains
    procedure :: InitTwoBodyOpChan
    procedure :: FinTwoBodyOpChan
    generic :: init => InitTwoBodyOpChan
    generic :: release => FinTwoBodyOpChan
  end type TwoBodyOpChan

  type :: TwoBodyOperator
    type(TwoBodyOpChan), allocatable :: MatCh(:,:)
    type(EleOrbits), pointer :: sps
    type(EleTwoBodySpace), pointer :: ms
    integer :: jrank = 0, prank = 1
  contains
    procedure :: InitTwoBodyOperator
    procedure :: FinTwoBodyOperator
    procedure :: SetTBME_scalar
    procedure :: SetTBME_tensor
    procedure :: GetTBME_scalar
    procedure :: GetTBME_tensor
    generic :: init => InitTwoBodyOperator
    generic :: fin => FinTwoBodyOperator
    generic :: SetTBME => SetTBME_scalar, SetTBME_tensor
    generic :: GetTBME => GetTBME_scalar, GetTBME_tensor
  end type TwoBodyOperator
contains

  subroutine InitTwoBodyOpChan( this, sps, chbra, chket )
    class(TwoBodyOpChan), intent(inout) :: this
    type(EleOrbits), intent(in), target :: sps
    type(EleTwoBodyChan), intent(in), target :: chbra, chket
    this%sps => sps
    this%chbra => chbra
    this%chket => chket
    call this%zeros( chbra%n_state, chket%n_state )
    this%zero = .false.
  end subroutine InitTwoBodyOpChan

  subroutine FinTwoBodyOpChan( this )
    class(TwoBodyOpChan), intent(inout) :: this
    if( this%zero ) return
    call this%fin()
    this%sps => null()
    this%chbra => null()
    this%chket => null()
  end subroutine FinTwoBodyOpChan

  subroutine InitTwoBodyOperator( this, ms, jrank, prank )
    class(TwoBodyOperator), intent(inout) :: this
    type(EleTwoBodySpace), intent(in), target :: ms
    integer, intent(in), optional :: jrank, prank
    integer :: ichbra, ichket
    type(EleTwoBodyChan), pointer :: chbra, chket
    integer :: jbra, jket, pbra, pket
    this%ms => ms
    this%sps => ms%sps
    if(present( jrank )) this%jrank = jrank
    if(present( prank )) this%prank = prank
    allocate( this%MatCh( ms%Nchan, ms%NChan) )
    do ichbra = 1, ms%NChan
      chbra => ms%jp(ichbra)
      jbra = chbra%j
      pbra = chbra%p
      do ichket = 1, ichbra
        chket => ms%jp(ichket)
        jket = chket%j
        pket = chket%p
        if( triag( jbra, jket, this%jrank )) cycle
        if( pbra * pket * this%prank == -1) cycle
        call this%MatCh(ichbra,ichket)%init( ms%sps, chbra, chket )
      end do
    end do
  end subroutine InitTwoBodyOperator

  subroutine FinTwoBodyOperator( this )
    class(TwoBodyOperator), intent(inout) :: this
    integer :: ichbra, ichket
    do ichbra = 1, this%ms%NChan
      do ichket = 1, ichbra
        call this%MatCh(ichbra,ichket)%release()
      end do
    end do
    this%ms => null()
    this%sps => null()
    this%jrank = 0
    this%prank = 1
  end subroutine FinTwoBodyOperator

  subroutine SetTBME_scalar( this, a, b, c, d, J, me )
    class(TwoBodyOperator), intent(inout) :: this
    integer, intent(in) :: a, b, c, d, J
    real(8), intent(in) :: me
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: ch, bra, ket, phase
    oa => this%sps%orb(a)
    ob => this%sps%orb(b)
    oc => this%sps%orb(c)
    od => this%sps%orb(d)
    if( mod(oa%l+ob%l,2) /= mod(oc%l+od%l,2) ) return
    if( triag( oa%j, ob%j, 2*J )) return
    if( triag( oc%j, od%j, 2*J )) return
    ch = this%ms%jp2idx(J,(-1)**(oa%l+ob%l))
    if( ch == 0 ) return
    bra = this%ms%jp(ch)%labels2n(a,b)
    ket = this%ms%jp(ch)%labels2n(c,d)
    phase = this%ms%jp(ch)%iphase(a,b) * this%ms%jp(ch)%iphase(c,d)
    this%MatCh(ch,ch)%m(bra,ket) = me * dble(phase)
    this%MatCh(ch,ch)%m(ket,bra) = me * dble(phase)
  end subroutine SetTBME_scalar

  subroutine SetTBME_tensor( this, a, b, c, d, Jab, Jcd, me )
    class(TwoBodyOperator), intent(inout) :: this
    integer, intent(in) :: a, b, c, d, Jab, Jcd
    real(8), intent(in) :: me
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: ichbra, ichket, bra, ket, phase, ibra, iket, iibra, iiket
    oa => this%sps%orb(a)
    ob => this%sps%orb(b)
    oc => this%sps%orb(c)
    od => this%sps%orb(d)
    if( (-1)**(oa%l+ob%l+oc%l+od%l)*this%prank == -1 ) return
    if( triag( oa%j, ob%j, 2*Jab )) return
    if( triag( oc%j, od%j, 2*Jcd )) return
    ichbra = this%ms%jp2idx(Jab,(-1)**(oa%l+ob%l))
    ichket = this%ms%jp2idx(Jcd,(-1)**(oc%l+od%l))
    if( ichbra * ichket == 0 ) return
    bra = this%ms%jp(ichbra)%labels2n(a,b)
    ket = this%ms%jp(ichket)%labels2n(c,d)
    if(bra * ket == 0) return
    phase = this%ms%jp(ichbra)%iphase(a,b) * this%ms%jp(ichket)%iphase(c,d)
    ibra = ichbra
    iket = ichket
    iibra = bra
    iiket = ket
    if( ichket > ichbra ) then
      ibra = ichket
      iket = ichbra
      iibra = ket
      iiket = bra
      phase = phase * (-1)**(this%jrank+Jab-Jcd)
    end if
    this%MatCh(ibra,iket)%m(iibra,iiket) = me * dble(phase)
    if( ibra /= iket ) return
    this%MatCh(ibra,iket)%m(iiket,iibra) = me * dble(phase)
  end subroutine SetTBME_tensor

  function GetTBME_scalar( this, a, b, c, d, J ) result(me)
    class(TwoBodyOperator), intent(in) :: this
    integer, intent(in) :: a, b, c, d, J
    real(8) :: me
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: ch, bra, ket, phase
    me = 0.d0
    oa => this%sps%orb(a)
    ob => this%sps%orb(b)
    oc => this%sps%orb(c)
    od => this%sps%orb(d)
    if( mod(oa%l+ob%l,2) /= mod(oc%l+od%l,2) ) return
    if( triag( oa%j, ob%j, 2*J )) return
    if( triag( oc%j, od%j, 2*J )) return
    ch = this%ms%jp2idx(J,(-1)**(oa%l+ob%l))
    if( ch == 0 ) return
    bra = this%ms%jp(ch)%labels2n(a,b)
    ket = this%ms%jp(ch)%labels2n(c,d)
    if( bra * ket == 0 ) return
    phase = this%ms%jp(ch)%iphase(a,b) * this%ms%jp(ch)%iphase(c,d)
    me = this%MatCh(ch,ch)%m(bra,ket) * dble(phase)
  end function GetTBME_scalar

  function GetTBME_tensor( this, a, b, c, d, Jab, Jcd ) result(me)
    class(TwoBodyOperator), intent(in) :: this
    integer, intent(in) :: a, b, c, d, Jab, Jcd
    real(8) :: me
    type(EleSingleParticleOrbit), pointer :: oa, ob, oc, od
    integer :: ichbra, ichket, bra, ket, phase, ibra, iket, iibra, iiket
    me = 0.d0
    oa => this%sps%orb(a)
    ob => this%sps%orb(b)
    oc => this%sps%orb(c)
    od => this%sps%orb(d)
    if( (-1)**(oa%l+ob%l+oc%l+od%l)*this%prank == -1 ) return
    if( triag( oa%j, ob%j, 2*Jab )) return
    if( triag( oc%j, od%j, 2*Jcd )) return
    ichbra = this%ms%jp2idx(Jab,(-1)**(oa%l+ob%l))
    ichket = this%ms%jp2idx(Jcd,(-1)**(oc%l+od%l))
    if( ichbra * ichket == 0 ) return
    bra = this%ms%jp(ichbra)%labels2n(a,b)
    ket = this%ms%jp(ichket)%labels2n(c,d)
    if(bra * ket == 0) return
    phase = this%ms%jp(ichbra)%iphase(a,b) * this%ms%jp(ichket)%iphase(c,d)
    ibra = ichbra
    iket = ichket
    iibra = bra
    iiket = ket
    if( ichket > ichbra ) then
      ibra = ichket
      iket = ichbra
      iibra = ket
      iiket = bra
      phase = phase * (-1)**(this%jrank+Jab-Jcd)
    end if
    me = this%MatCh(ibra,iket)%m(iibra,iiket) * dble(phase)
  end function GetTBME_tensor

end module TwoBodyTerms
