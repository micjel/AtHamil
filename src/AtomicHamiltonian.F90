module AtomicHamiltonian
  use omp_lib
  use LinAlgLib
  use ElectronTwoBodySpace
  use EleSingleParticleState
  use OneBodyTerms, only: AtomicHamilChan
  implicit none

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

  real(8), private :: a_0
contains
  subroutine FinAtomicHamil(this)
    class(AtomicHamil), intent(inout) :: this
    integer :: chbra, chket
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

  subroutine SetAtomicHamil(this, NMesh, rmax, NMesh_mom, pmax)
    use OneBodyTerms
    use TwoBodyCoulomb
    use TwoBodyDarwin
    use TwoBodySpinContact
    use TwoBodyOrbitOrbit
    use TwoBodySpinOrbit
    use TwoBodySpinDipole
    class(AtomicHamil), intent(inout) :: this
    integer, intent(in) :: NMesh, NMesh_Mom
    real(8), intent(in) :: rmax, pmax
    real(8) :: ti

    ti = omp_get_wtime()
    select case(this%ms%sps%basis)
    case("HO", "ho")
      call set_kinetic_ho( this%kinetic, this%ms%sps, this%ms%zeta )
      call set_coulomb_ho( this%potential, this%ms%sps, this%ms%zeta, NMesh, rmax )
    case("AO", "ao")
      call set_kinetic_hydrogen( this%kinetic, this%ms%sps, this%ms%zeta, NMesh_mom, pmax )
      call set_coulomb_hydrogen( this%potential, this%ms%sps, this%ms%zeta, NMesh, rmax )
    case("LO", "lo", "LO-2nl", "lo-2nl")
      call set_kinetic_laguerre( this%kinetic, this%ms%sps, this%ms%zeta )
      call set_coulomb_laguerre( this%potential, this%ms%sps, this%ms%zeta, NMesh, rmax )
      call set_kinetic_correction_laguerre( this%kinetic_p4, this%ms%sps, this%ms%zeta, NMesh_mom, pmax )
      call set_Darwin_laguerre( this%Darwin_term, this%ms%sps, this%ms%zeta )
      call set_spin_orbit_laguerre( this%LS_term, this%ms%sps, this%ms%zeta, NMesh, rmax )

      call set_ee_coulomb_laguerre( this%ee_coulomb, this%ms, NMesh, rmax )
      call set_ee_darwin_laguerre( this%ee_darwin, this%ms, NMesh, rmax )
      call set_ee_spin_contact_laguerre( this%ee_spin_contact, this%ms, NMesh, rmax )
      call set_ee_orbit_orbit_laguerre( this%ee_orbit_orbit, this%ms, NMesh, rmax )
      call set_ee_spin_orbit_laguerre( this%ee_spin_orbit, this%ms, NMesh, rmax )
      call set_ee_spin_dipole_laguerre( this%ee_spin_dipole, this%ms, NMesh, rmax )
    end select
    write(*,"(a,f12.6,a)") "SetHamil:  ", omp_get_wtime() - ti, " sec"
  end subroutine SetAtomicHamil

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

end module AtomicHamiltonian
