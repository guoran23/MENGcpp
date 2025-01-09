!2023/03/31 from spline-gkx zhixin.lu@ipp.mpg.de
module equil_class
  use mpi_class
  use bspline_module
  use util_math

  implicit none
!#IFNDEF rank
! integer::rank=0
!#ENDIF
  type equil_cls
    private
!--1.0 controls
    integer iequmodel
    logical set_zerof  !--set fpol=0 manually
    real*8 c1adhoc,c2adhoc,rmaxis_adhoc,Bmaxis_adhoc
!--1.1 eqdsk
    character*100 fname
    character*100 profilename
    character*10 case(6)
    integer idum, nw, nh, nbbbs, limitr
    !real*8 rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum
    real*8 zdim, rcentr, rleft, zmid, zmaxis, simag, sibry, bcentr, current, xdum
    real*8,public::rmaxis,rdim
    real*8, allocatable, dimension(:)::fpol, pres, ffprim, pprime, qpsi, rbbbs, zbbbs, rlim, zlim
    real*8, allocatable, dimension(:, :)::psirz
!--1.2 derived variables
!--derived variables
    real*8,public::Bmaxis ! added 20200825. B on axis. bcentr--rcentr; Bmaxis--rmaxis
    real*8,public::Bref = 1
    real*8 signBphi
    real*8 dr, dz, dpsi, psiwid
    real*8, allocatable, dimension(:)::r1d, z1d, psi1d
    integer isetpsi
    real*8 fBphi
!--1.3 spline class
    logical extrap
    integer r_order_bsp, z_order_bsp, phi_order_bsp, psi_order_bsp
    integer r_order_bsp3d, z_order_bsp3d, phi_order_bsp3d
!--bspline coefficients in (R,Z)
    type(bspline_2d)::radrz_bsp

    type(bspline_1d)::F1d_bsp
    type(bspline_1d)::q1d_bsp
!--bspline in (rad,the)
    type(bspline_1d)::rbbbs_bsp, zbbbs_bsp
    type(bspline_2d)::rrt_bsp, zrt_bsp, Brt_bsp
    type(bspline_2d)::Brad_co_bsp, Bthe_co_bsp
    type(bspline_2d)::curlbrad_ct_bsp, curlbthe_ct_bsp, curlbphi_ct_bsp
    integer::nrad_sp = 32, nthe_sp = 64
    real*8::radmax_sp = 0.7d0, radmin_sp = 0.3d0
    real*8 drad_sp, dthe_sp
    real*8, dimension(:), allocatable::rad1d_sp, the1d_sp
!--bspline in (rad,the,phi)
    type(bspline_3d)::eta3d_bsp
    type(bspline_3d)::the3d_bsp
    real*8::nphieachfa = 4 !note:can be non-integer; set to 4 for cubic
    integer::nphi = 7 !toroidal grid #
    integer::nphiMult = 2 !torus wedge factor
    real*8::phimaxfull = twopi, phimin = 0.d0
    real*8 phimax, phiwid, dphi
    real*8 phimin_fa, phimax_fa, dphi_fa
    integer::nphifa_sp = 8 !toroidal grid # of 3d Spline
    integer::nintg_fa = 10
    real*8, dimension(:), allocatable::phifa1d_sp
!--3.1 physics variables
    !1:RhoN,BetaN,BRef->NREF,TREF;2:NREF,TREF,BRef->RhoN,BetaN;T:eV;n:m^3
    integer::isource_ref = 1
    real*8,public::rhoN, betaN, nref, Tref
!--for diagnosis
!--methods--------------
  contains
    procedure::equil_cls_init
    procedure::equil_cls_readinput
    procedure::equil_cls_calc_derived
    !(R,Z) functions
    procedure::getradRZ
    procedure::gettheRZ
    !(rad,the) functions
    procedure::getBdirect
    procedure::getBrad_co_direct
    procedure::getBthe_co_direct
    !procedure::getBrz
    !procedure::getdBdRrz
    !procedure::getdBdZrz
    procedure::getB
    procedure::getdBdrad
    procedure::getdBdthe
    procedure::getBthe_ct
    procedure::getBphi_ct
    procedure::getcurlbrad_ct
    procedure::getcurlbthe_ct
    procedure::getcurlbphi_ct

    !---for Rho* terms of GC motion---
    procedure::getBrad_co
    procedure::getBthe_co
    procedure::getBphi_co
    procedure::getdBrad_codrad
    procedure::getdBrad_codthe
    procedure::getdBthe_codrad
    procedure::getdBthe_codthe
    procedure::calc_curlb_ct_direct

    procedure::calcBJgij

    procedure::getgij_fa
    procedure::getgij_rtp

    procedure::getR
    procedure::getZ
    procedure::getdRdrad
    procedure::getdZdrad
    procedure::getdRdthe
    procedure::getdZdthe
    procedure::calcdrtdRZ
    procedure::getjaco2
    procedure::getjaco3

    procedure::getqloc_rt
    procedure::getdqdrloc_rt
    procedure::getfpol_r
    procedure::getdfpoldr_r
    procedure::getdpsidr
    procedure::getpsi

    procedure::geteta
    procedure::getdetadrad
    procedure::getdetadthe
    procedure::getdetadphi

    procedure::getthe3d
    procedure::getdthe3ddrad
    procedure::getdthe3ddthe
    procedure::getdthe3ddphi

    procedure::getfeq

    procedure::equil_cls_showinfo

  end type equil_cls
!----------------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------------
  subroutine equil_cls_init(this)
    implicit none
    class(equil_cls)::this
!--idx,idy: derivative order
    integer neqdsk, fic, fjc, fkc, iflag, idx, idy, idpsi, nfile,stat
    real*8 radval, theval, rval, zval, psival, fval1, fval2, fval3
    real*8, dimension(:), allocatable::thebbbs, thebbbs_q, rbbbs_q, zbbbs_q !_q:for interpolation
    real*8, dimension(:), allocatable::thebbbs_q_cp, rbbbs_q_cp, zbbbs_q_cp !_q:for interpolation
    real*8, dimension(:, :), allocatable::rval2d, zval2d, fval2d1_rz, fval2d1_rt, fval2d2_rt, fval2d3_rt
    logical extrap
    real*8 ytmp
    logical ifswap
    integer isvtmp, isvmax
    real*8 r0, r1, r2, z0, z1, z2, relerr
    integer nshift, nsame
!--eta3d--
    integer nintg_tot
    real*8 phiend
    real*8, dimension(:, :, :), allocatable::fval3d_rtp
    real*8, dimension(:, :, :), allocatable::fval3d_rep

    real*8, dimension(:), allocatable::xx1d, yy1d
    !--Curl of b--
    real*8 ftmp1, ftmp2, ftmp3
!--
    if (rank .eq. 0) write (*, *) '===========================equil_class init starts============================'
!--1. read eqdsk
    call this%equil_cls_readinput
    call this%equil_cls_calc_derived
    call this%equil_cls_showinfo
!--1.1
   if(this%iequmodel.eq.1)then
    neqdsk = 111
    if (rank .eq. 0) write (*, *) '----eqdsk file: ', trim(this%fname)
    open (neqdsk, file=trim(this%fname), status='old', READONLY,iostat=stat)
    if(stat.ne.0)write(*,*)' open eqdsk file fails '
    
    read (neqdsk, 2000) (this%case(fic), fic=1, 6), this%idum, this%nw, this%nh
    read (neqdsk, 2020) this%rdim, this%zdim, this%rcentr, this%rleft, this%zmid
    read (neqdsk, 2020) this%rmaxis, this%zmaxis, this%simag, this%sibry, this%bcentr
    read (neqdsk, 2020) this%current, this%simag, this%xdum, this%rmaxis, this%xdum
    read (neqdsk, 2020) this%zmaxis, this%xdum, this%sibry, this%xdum, this%xdum

!--note psirz(nw,nh), i.e., (r,z), in G_EQDSK, L.Lao !!!
    allocate (this%fpol(this%nw), this%pres(this%nw), this%ffprim(this%nw), this%pprime(this%nw), this%qpsi(this%nw))
    allocate (this%psirz(this%nw, this%nh))

    read (neqdsk, 2020) (this%fpol(fic), fic=1, this%nw)
    read (neqdsk, 2020) (this%pres(fic), fic=1, this%nw)
    read (neqdsk, 2020) (this%ffprim(fic), fic=1, this%nw)
    read (neqdsk, 2020) (this%pprime(fic), fic=1, this%nw)
    read (neqdsk, 2020) ((this%psirz(fic, fjc), fic=1, this%nw), fjc=1, this%nh)
    read (neqdsk, 2020) (this%qpsi(fic), fic=1, this%nw)
    read (neqdsk, 2022) this%nbbbs, this%limitr

    if (rank .eq. 0) write (*, *) '----nbbbs,limitr: ', this%nbbbs, this%limitr
    allocate (this%rbbbs(this%nbbbs), this%zbbbs(this%nbbbs), this%rlim(this%limitr), this%zlim(this%limitr))

    read (neqdsk, 2020) (this%rbbbs(fic), this%zbbbs(fic), fic=1, this%nbbbs)
    read (neqdsk, 2020) (this%rlim(fic), this%zlim(fic), fic=1, this%limitr)

    close (neqdsk)
!--manual operation
    if (this%set_zerof) then
      if (rank .eq. 0) write (*, *) '-------- fpol set to constant (ffprim=0) manually! --------'
      do fic = 1, this%nw
        this%fpol(fic) = this%fpol(1)
        this%ffprim(fic) = 0.d0
      end do
    end if
    if (this%fBphi .ne. 1.d0) then
      if (rank .eq. 0) write (*, '(A50,e15.5)') '-------- fpol multiplied by fBphi=', this%fBphi
      this%fpol(:) = this%fpol(:)*this%fBphi
      this%ffprim(:) = this%ffprim(:)*this%fBphi
    end if
    if (this%isetpsi .eq. 1) then
      !make axis psi to zero
      if (rank .eq. 0) write (*, *) ' shift psirz,simag,sibry by -simag '
      this%psirz(:, :) = this%psirz(:, :) - this%simag
      this%sibry = this%sibry - this%simag
      this%simag = 0.d0
      !make simag<sibry
      if (this%sibry .lt. this%simag) then
        if (rank .eq. 0) write (*, *) ' invert psirz,simag,sibry, ffprim '
        this%sibry = -this%sibry
        this%simag = -this%simag
        this%psirz(:, :) = -this%psirz(:, :)
        this%ffprim(:) = -this%ffprim(:)
      end if
      if (rank .eq. 0) write (*, '(A20,e12.5,e12.5)') ' simag,sibry=', this%simag, this%sibry
    end if
!--derived variables
    this%Bmaxis = this%fpol(1)/this%rmaxis
!--check data
    if (rank .eq. 0) then
      !write(*,*)'----eqdsk data:'
      !write(*,2000) (this%case(fic),fic=1,6),this%idum,this%nw,this%nh
      !write(*,2020) this%rdim,this%zdim,this%rcentr,this%rleft,this%zmid
      !write(*,2020) this%rmaxis,this%zmaxis,this%simag,this%sibry,this%bcentr
      !write(*,2020) this%current,this%simag,this%xdum,this%rmaxis,this%xdum
      !write(*,2020) this%zmaxis,this%xdum,this%sibry,this%xdum,this%xdum
!   write(*,2020) this%fpol
!   write(*,2020) (this%rbbbs(fic),this%zbbbs(fic),fic=1,this%nbbbs)
      write (*, *) '...... end of eqdsk data------'
    end if
2000 format(6a8, 3i4)
2020 format(5e16.9)
2022 format(2i5)

!--2. calc derived variables
    this%signBphi = dsign(1.d0, this%fpol(1))
    allocate (this%r1d(this%nw), this%z1d(this%nh), this%psi1d(this%nw))
!--check the definition of r1d,z1d
    this%dr = this%rdim/(this%nw - 1)
    this%dz = this%zdim/(this%nh - 1)
    this%dpsi = (this%sibry - this%simag)/(this%nw - 1)
    this%psiwid = this%sibry - this%simag

    do fic = 1, this%nw
      this%r1d(fic) = this%rleft + this%dr*(fic - 1)
    end do
    do fic = 1, this%nh
      this%z1d(fic) = this%zmid - this%zdim/2.0 + this%dz*(fic - 1)
    end do
    do fic = 1, this%nw
      this%psi1d(fic) = this%simag + this%dpsi*(fic - 1)
    end do
    if (rank .eq. 0) then
  write(*,*)'----rmaxis,rleft,rdim,zmid,zdim,psimag,psibry:',this%rmaxis,this%rleft,this%rdim,this%zmid,this%zdim,this%simag,this%sibry
      write (*, *) ' Bmaxis=', this%Bmaxis, ', zmaxis=', this%zmaxis
!  write(*,*)'----F(0)=',this%fpol(1)
!  write(*,*)'----dr,dz,dpsi:',this%dr,this%dz,this%dpsi
!  write(*,*)'----r1d:',this%r1d(1:3),' ... ',this%r1d(this%nw-2:this%nw)
!  write(*,*)'----z1d:',this%z1d(1:3),' ... ',this%z1d(this%nh-2:this%nh)
      write (*, *) '----psi1d:', this%psi1d(1:3), ' ... ', this%psi1d(this%nw - 2:this%nw)
    end if
!--prepare spline: moved to readinput
    extrap = this%extrap
!  this%r_order_bsp=4
!  this%z_order_bsp=4
!  this%phi_order_bsp=4
!  this%psi_order_bsp=4
!
!  this%r_order_bsp3d=3
!  this%z_order_bsp3d=3
!  this%phi_order_bsp3d=3
    if (rank .eq. 0) then
    write(*,'(A20,4I5,A20,3I5,A10,L6)')'r,z,phi,psi_order_bsp=',this%r_order_bsp,this%z_order_bsp,this%phi_order_bsp,this%psi_order_bsp, &
        ', r,z,phi_order_bsp3d=', this%r_order_bsp3d, this%z_order_bsp3d, this%phi_order_bsp3d, ',extrap=', this%extrap
    end if
!----calc. BR,BZ,B
    allocate (rval2d(this%nw, this%nh), zval2d(this%nw, this%nh), fval2d1_rz(this%nw, this%nh))
!--1. RADRZ
    fval2d1_rz(:, :) = sqrt((this%psirz(:, :) - this%simag)/(this%sibry - this%simag))
    call this%radrz_bsp%initialize(this%r1d, this%z1d, fval2d1_rz, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in radrz bsp:', iflag
!--calc. theta
!  idx=0
!  idy=1
!  do fjc=1,this%nh
!   do fic=1,this%nw
!    rval=this%r1d(fic)
!    zval=this%z1d(fjc)
!    fval2d_rz(fic,fjc)=datan2(zval-this%zmaxis,rval-this%rmaxis)
!   enddo
!  enddo
!  call this%theta2d_bsp%initialize(this%r1d,this%z1d,fval2d_rz,this%r_order_bsp,this%z_order_bsp,iflag,extrap)
!--5.1. 1d F interpolation for calc. Bphi
!--note consistence with inpoly!!!
    if (this%simag .lt. this%sibry) then
      call this%F1d_bsp%initialize(this%psi1d, this%fpol, this%psi_order_bsp, iflag, extrap)
      call this%q1d_bsp%initialize(this%psi1d, this%qpsi, this%psi_order_bsp, iflag, extrap)
    else
      allocate (xx1d(this%nw), yy1d(this%nw))
      do fic = 1, this%nw
        xx1d(this%nw + 1 - fic) = this%psi1d(fic)
        yy1d(this%nw + 1 - fic) = this%fpol(fic)
      end do
      call this%F1d_bsp%initialize(xx1d, yy1d, this%psi_order_bsp, iflag, extrap)
      do fic = 1, this%nw
        xx1d(this%nw + 1 - fic) = this%psi1d(fic)
        yy1d(this%nw + 1 - fic) = this%qpsi(fic)
      end do
      call this%q1d_bsp%initialize(xx1d, yy1d, this%psi_order_bsp, iflag, extrap)
    end if
!--B
    !do fjc=1,this%nh
    ! do fic=1,this%nw
    !  rval=this%r1d(fic)
    !  zval=this%z1d(fjc)
    !  idx=1
    !  idy=0
    !  call this%radrz_bsp%evaluate(r1,z1,idx,idy,fval1,iflag)
    !  idx=0
    !  idy=1
    !  call this%radrz_bsp%evaluate(r1,z1,idx,idy,fval2,iflag)
    !  fval3=fval2d1_rz(fic,fjc)
    !  fval2d2_rz(fic,fjc)=sqrt( this%getfpol_r(fval3)**2 &
    !          +(fval1**2+fval2**2)*this%getdpsidr(fval3)**2 )/rval
    ! enddo
    !enddo

    !call this%Brz_bsp%initialize(this%r1d,this%z1d,fval2d2_rz,this%r_order_bsp,this%z_order_bsp,iflag,extrap)
    !if(rank.eq.0.and.iflag.ne.0)write(*,*)'****Error in Brz bsp:',iflag
!--output
    if (rank == 0) then
      nfile = 201
      open (nfile, file='rzbbbs.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open rzbbbs.txt fails '
      do fic = 1, this%nbbbs
        write (nfile, '(2e16.4)') this%rbbbs(fic), this%zbbbs(fic)
      end do
      close (nfile)

      nfile = 201
      open (nfile, file='rzlim.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open rzlim.txt fails '
      do fic = 1, this%limitr
        write (nfile, '(2e16.4)') this%rlim(fic), this%zlim(fic)
      end do
      close (nfile)
    end if
    if (rank == 0) then
      nfile = 201
      open (nfile, file='radrz.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open radrz.txt fails '
      write (nfile, '(1e16.9)') fval2d1_rz
      close (nfile)
    end if
!  if(rank==0) then
!   open(nfile,file='Brz.txt')
!   write(nfile,'(1e16.9)')fval2d2_rz
!   close(nfile)
!  endif

    !calc. R(rad,the),Z(rad,the)
  allocate(thebbbs(this%nbbbs))
  allocate(this%rad1d_sp(this%nrad_sp),this%the1d_sp(this%nthe_sp))
  allocate(fval2d1_rt(this%nrad_sp,this%nthe_sp),fval2d2_rt(this%nrad_sp,this%nthe_sp),fval2d3_rt(this%nrad_sp,this%nthe_sp))

    thebbbs(:) = 0.d0

    do fic = 1, this%nbbbs
      thebbbs(fic) = atan2(this%zbbbs(fic) - this%zmaxis, this%rbbbs(fic) - this%rmaxis)
      if (thebbbs(fic) .lt. 0.d0) thebbbs(fic) = thebbbs(fic) + twopi
    end do
    !sort theta and r,z along boundary
    ifswap = .true.
    do while (ifswap)
      ifswap = .false.
      do fic = 1, this%nbbbs - 1
        if (thebbbs(fic) .gt. thebbbs(fic + 1)) then
          ytmp = thebbbs(fic)
          thebbbs(fic) = thebbbs(fic + 1)
          thebbbs(fic + 1) = ytmp
          ytmp = this%rbbbs(fic)
          this%rbbbs(fic) = this%rbbbs(fic + 1)
          this%rbbbs(fic + 1) = ytmp
          ytmp = this%zbbbs(fic)
          this%zbbbs(fic) = this%zbbbs(fic + 1)
          this%zbbbs(fic + 1) = ytmp
          ifswap = .true.
        end if
      end do
    end do

    allocate (thebbbs_q(this%nbbbs + 2), rbbbs_q(this%nbbbs + 2), zbbbs_q(this%nbbbs + 2))
    allocate (thebbbs_q_cp(this%nbbbs + 2), rbbbs_q_cp(this%nbbbs + 2), zbbbs_q_cp(this%nbbbs + 2))
    thebbbs_q(:) = 0.d0
    rbbbs_q(:) = 0.d0
    zbbbs_q(:) = 0.d0

    thebbbs_q(2:this%nbbbs + 1) = thebbbs
    if (thebbbs_q(this%nbbbs + 1) - twopi .eq. thebbbs_q(2)) then
      nshift = 1
    else
      nshift = 0
    end if
    !add redudant point in head and end for cover theta=0,2pi
    thebbbs_q(1) = thebbbs_q(this%nbbbs + 1 - nshift) - twopi
    thebbbs_q(this%nbbbs + 2) = thebbbs_q(2 + nshift) + twopi

    rbbbs_q(2:this%nbbbs + 1) = this%rbbbs
    rbbbs_q(1) = rbbbs_q(this%nbbbs + 1 - nshift)
    rbbbs_q(this%nbbbs + 2) = rbbbs_q(2 + nshift)

    zbbbs_q(2:this%nbbbs + 1) = this%zbbbs
    zbbbs_q(1) = zbbbs_q(this%nbbbs + 1 - nshift)
    zbbbs_q(this%nbbbs + 2) = zbbbs_q(2 + nshift)
    !avoid same theta and error in bsp init
    nsame = 0
    thebbbs_q_cp(1) = thebbbs_q(1)
    rbbbs_q_cp(1) = rbbbs_q(1)
    zbbbs_q_cp(1) = zbbbs_q(1)
    do fic = 2, this%nbbbs + 2
      if (thebbbs_q(fic) == thebbbs_q(fic - 1)) then
        nsame = nsame + 1
        if (rank .eq. 0) write (*, *) 'same value at ', fic - 1
      else
        thebbbs_q_cp(fic - nsame) = thebbbs_q(fic)
        rbbbs_q_cp(fic - nsame) = rbbbs_q(fic)
        zbbbs_q_cp(fic - nsame) = zbbbs_q(fic)
      end if
    end do
    if (nsame .ge. 1) then
      if (rank .eq. 0) write (*, *) 'same theta-(R,Z), nsame=', nsame
      deallocate (thebbbs_q, rbbbs_q, zbbbs_q)
      allocate (thebbbs_q(this%nbbbs + 2 - nsame), rbbbs_q(this%nbbbs + 2 - nsame), zbbbs_q(this%nbbbs + 2 - nsame))
      thebbbs_q(:) = thebbbs_q_cp(1:this%nbbbs + 2 - nsame)
      rbbbs_q(:) = rbbbs_q_cp(1:this%nbbbs + 2 - nsame)
      zbbbs_q(:) = zbbbs_q_cp(1:this%nbbbs + 2 - nsame)
    end if
    !bsp init
    call this%rbbbs_bsp%initialize(thebbbs_q, rbbbs_q, this%psi_order_bsp, iflag, extrap)
    if (iflag .ne. 0 .and. rank .eq. 0) write (*, *) 'Error in rbbbs_bsp init:', iflag
    call this%zbbbs_bsp%initialize(thebbbs_q, zbbbs_q, this%psi_order_bsp, iflag, extrap)
    if (iflag .ne. 0 .and. rank .eq. 0) write (*, *) 'Error in zbbbs_bsp init:', iflag

    this%drad_sp = (this%radmax_sp - this%radmin_sp)/(this%nrad_sp - 1)
    this%dthe_sp = twopi/(this%nthe_sp - 1)
    do fic = 1, this%nrad_sp
      this%rad1d_sp(fic) = this%radmin_sp + (fic - 1)*this%drad_sp
    end do
    do fic = 1, this%nthe_sp
      this%the1d_sp(fic) = (fic - 1)*this%dthe_sp
    end do

    fval2d1_rt = 0.d0
    fval2d2_rt = 0.d0
    fval2d3_rt = 0.d0

    idx = 0
    idy = 0
    do fic = 1, this%nthe_sp
      theval = this%the1d_sp(fic)
      call this%rbbbs_bsp%evaluate(theval, idx, fval1, iflag)
      fval2d1_rt(this%nrad_sp, fic) = fval1

      call this%zbbbs_bsp%evaluate(theval, idx, fval1, iflag)
      fval2d2_rt(this%nrad_sp, fic) = fval1
      !force periodicity at 2pi
!   if(fic.eq.this%nthe_sp)fval2d1_rt(this%nrad_sp,fic)=fval2d1_rt(this%nrad_sp,1)
!   if(fic.eq.this%nthe_sp)fval2d2_rt(this%nrad_sp,fic)=fval2d2_rt(this%nrad_sp,1)

      relerr = 1d-5/sqrt(1.d0*this%nthe_sp**2 + this%nrad_sp**2)
      isvmax = 100
      do fjc = 1, this%nrad_sp - 0
        if (fjc .ge. 2) then
          r0 = fval2d1_rt(fjc - 1, fic)
          z0 = fval2d2_rt(fjc - 1, fic)
        elseif (fjc .eq. 1) then
          r0 = this%rmaxis
          z0 = this%zmaxis
        end if
        r2 = fval2d1_rt(this%nrad_sp, fic)
        z2 = fval2d2_rt(this%nrad_sp, fic)

        isvtmp = 0
        fval1 = 1.d0
        do while (abs(fval1) .gt. relerr .and. isvtmp .le. isvmax)
          isvtmp = isvtmp + 1
          r1 = (r0 + r2)/2
          z1 = (z0 + z2)/2
          call this%radrz_bsp%evaluate(r1, z1, idx, idy, fval1, iflag)
          !if(iflag.ne.0)write(*,*)'Error in radrz_bsp:',iflag

          fval1 = fval1 - this%rad1d_sp(fjc)
          if (fval1 .lt. 0) then
            r0 = r1
            z0 = z1
          else
            r2 = r1
            z2 = z1
          end if
        end do !WHILE
        fval2d1_rt(fjc, fic) = r1
        fval2d2_rt(fjc, fic) = z1
      end do !fjc-rad
!  write(*,*)'do #:',isvtmp
    end do !fic-the

    call this%rrt_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d1_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:rrt bsp init:', iflag

    call this%zrt_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d2_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:zrt bsp init:', iflag

    !--
    do fic = 1, this%nrad_sp
      do fjc = 1, this%nthe_sp
        !need F,drzdRZ,R in B=sqrt( (this%getdpsidr(rad)**2*g11+FF**2)/RR**2 )
        fval2d3_rt(fic, fjc) = this%getBdirect(this%rad1d_sp(fic), this%the1d_sp(fjc))
      end do
    end do
    call this%Brt_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d3_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:Brt bsp init:', iflag

    !--construct fa3d--
    this%phiwid = this%phimaxfull - this%phimin
    this%phimax = this%phiwid/this%nphiMult
    this%dphi = this%phiwid/this%nphi
    this%phimin_fa = -this%dphi*this%nphieachfa/2
    this%phimax_fa = this%dphi*this%nphieachfa/2
    this%dphi_fa = this%dphi*this%nphieachfa/(this%nphifa_sp - 1)

    allocate (this%phifa1d_sp(this%nphifa_sp))
    allocate (fval3d_rtp(this%nrad_sp, this%nthe_sp, this%nphifa_sp))
    allocate (fval3d_rep(this%nrad_sp, this%nthe_sp, this%nphifa_sp))

    do fic = 1, this%nphifa_sp
      this%phifa1d_sp(fic) = this%phimin_fa + (fic - 1)*this%dphi_fa
    end do

    do fic = 1, this%nrad_sp
      do fjc = 1, this%nthe_sp
        do fkc = 1, this%nphifa_sp
          nintg_tot = this%nintg_fa*(1 + floor(abs(this%phifa1d_sp(fkc)/this%dphi_fa)))
          phiend = -this%phifa1d_sp(fkc)
          fval3d_rtp(fic, fjc, fkc) = onestep_dtheodphi_nstep(this, this%rad1d_sp(fic), this%the1d_sp(fjc), phiend, nintg_tot)
          phiend = +this%phifa1d_sp(fkc)
          fval3d_rep(fic, fjc, fkc) = onestep_dtheodphi_nstep(this, this%rad1d_sp(fic), this%the1d_sp(fjc), phiend, nintg_tot)
        end do
      end do
    end do

    call this%eta3d_bsp%initialize(this%rad1d_sp, this%the1d_sp, this%phifa1d_sp, fval3d_rtp, &
                                   this%r_order_bsp3d, this%z_order_bsp3d, this%phi_order_bsp3d, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:eta3d bsp init:', iflag
    call this%the3d_bsp%initialize(this%rad1d_sp, this%the1d_sp, this%phifa1d_sp, fval3d_rep, &
                                   this%r_order_bsp3d, this%z_order_bsp3d, this%phi_order_bsp3d, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:the3d bsp init:', iflag

    !--Visualize for debug--
    if (rank == 0) then
      write(*,*)' write rzbbbs_sort,rzbbbs_sp,RZ2dsp.txt'
      nfile = 201
      open (nfile, file='rzbbbs_sort.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open rzbbbs_sort.txt fails '
      do fic = 1, this%nbbbs
        write (nfile, '(2e16.4)') this%rbbbs(fic), this%zbbbs(fic)
      end do
      close (nfile)

      nfile = 201
      open (nfile, file='rzbbbs_sp.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open rzbbbs_sp.txt fails '
      do fic = 1, size(thebbbs_q)
        write (nfile, '(2e16.4)') rbbbs_q(fic), zbbbs_q(fic), thebbbs_q(fic)
      end do
      close (nfile)

      nfile = 201
      open (nfile, file='RZ2dsp.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open RZ2dsp.txt fails '
      do fjc = 1, this%nthe_sp
        do fic = 1, this%nrad_sp
          write (nfile, '(20e16.4)') fval2d1_rt(fic, fjc), fval2d2_rt(fic, fjc), fval2d3_rt(fic, fjc)
        end do
      end do
      close (nfile)

      nfile = 201
      open (nfile, file='eta3dsp.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open eta3dsp.txt fails '
      do fkc = 1, this%nphifa_sp
        do fjc = 1, this%nthe_sp
          do fic = 1, this%nrad_sp
            write (nfile, '(2e16.4)') fval3d_rtp(fic, fjc, fkc)
          end do
        end do
      end do
      close (nfile)

      nfile = 201
      open (nfile, file='the3dsp.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open the3dsp.txt fails '
      do fkc = 1, this%nphifa_sp
        do fjc = 1, this%nthe_sp
          do fic = 1, this%nrad_sp
            write (nfile, '(2e16.4)') fval3d_rep(fic, fjc, fkc)
          end do
        end do
      end do
      close (nfile)
    end if

    !--co-variant B components
    if (rank .eq. 0) write (*, *) '--------calc. co-variant B components--------'
    do fic = 1, this%nrad_sp
      do fjc = 1, this%nthe_sp
        fval2d1_rt(fic, fjc) = this%getBrad_co_direct(this%rad1d_sp(fic), this%the1d_sp(fjc))
        fval2d2_rt(fic, fjc) = this%getBthe_co_direct(this%rad1d_sp(fic), this%the1d_sp(fjc))
      end do
    end do

    call this%Brad_co_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d1_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:Brad_co bsp init:', iflag

    call this%Bthe_co_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d2_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:Bthe_co bsp init:', iflag

    !--Visualize Brt_co for debug--
    if (rank == 0) then
      nfile = 201
      open (nfile, file='Brt_co2d.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open Brt_co2d.txt fails '
      do fjc = 1, this%nthe_sp
        do fic = 1, this%nrad_sp
          write (nfile, '(20e16.4)') fval2d1_rt(fic, fjc), fval2d2_rt(fic, fjc)
        end do
      end do
      close (nfile)
    end if

    !--Calculate curl of b--
    if (rank .eq. 0) write (*, *) '--------calc. curl{b}--------'
    do fic = 1, this%nrad_sp
      do fjc = 1, this%nthe_sp
        call this%calc_curlb_ct_direct(this%rad1d_sp(fic), this%the1d_sp(fjc), ftmp1, ftmp2, ftmp3)
        fval2d1_rt(fic, fjc) = ftmp1
        fval2d2_rt(fic, fjc) = ftmp2
        fval2d3_rt(fic, fjc) = ftmp3
      end do
    end do
    !--Set Bsplines--
   call this%curlbrad_ct_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d1_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:curlbrad_ct bsp init:', iflag

   call this%curlbthe_ct_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d2_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:curlbthe_ct bsp init:', iflag

   call this%curlbphi_ct_bsp%initialize(this%rad1d_sp, this%the1d_sp, fval2d3_rt, this%r_order_bsp, this%z_order_bsp, iflag, extrap)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) 'Error:curlbphi_ct bsp init:', iflag

    !--Visualize curl{b} for debug--
    if (rank == 0) then
      nfile = 201
      open (nfile, file='curlb_ct.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open curlb_ct.txt fails '
      do fjc = 1, this%nthe_sp
        do fic = 1, this%nrad_sp
          write (nfile, '(20e16.4)') fval2d1_rt(fic, fjc), fval2d2_rt(fic, fjc), fval2d3_rt(fic, fjc)
        end do
      end do
      close (nfile)
    end if

   elseif(this%iequmodel==2)then
    call equil_cls_init_adhoc(this)
   endif
    !--Visualize--
    call record_var2d(this)
    call record_var3d(this)

    if (rank .eq. 0) write (*, *) '===========================equil_class init ends============================'

  end Subroutine equil_cls_init

!----------------------------------------------------------------------------------
  subroutine equil_cls_init_adhoc(this)
    implicit none
    class(equil_cls)::this
    integer fic

    !--Set variables--
    if(rank.eq.0)write(*,*)'----initialize iequmodel==2----'
    this%rmaxis=this%rmaxis_adhoc
    this%Bmaxis=this%Bmaxis_adhoc
    this%zmaxis=0.d0

    !--Set arrays for diagnosis equvar2d,3d
    allocate(this%rad1d_sp(this%nrad_sp),this%the1d_sp(this%nthe_sp))

    this%drad_sp = (this%radmax_sp - this%radmin_sp)/(this%nrad_sp - 1)
    this%dthe_sp = twopi/(this%nthe_sp - 1)
    do fic = 1, this%nrad_sp
      this%rad1d_sp(fic) = this%radmin_sp + (fic - 1)*this%drad_sp
    end do
    do fic = 1, this%nthe_sp
      this%the1d_sp(fic) = (fic - 1)*this%dthe_sp
    end do

    !--construct fa3d--
    this%phiwid = this%phimaxfull - this%phimin
    this%phimax = this%phiwid/this%nphiMult
    this%dphi = this%phiwid/this%nphi
    this%phimin_fa = -this%dphi*this%nphieachfa/2
    this%phimax_fa = this%dphi*this%nphieachfa/2
    this%dphi_fa = this%dphi*this%nphieachfa/(this%nphifa_sp - 1)

    allocate (this%phifa1d_sp(this%nphifa_sp))

    do fic = 1, this%nphifa_sp
      this%phifa1d_sp(fic) = this%phimin_fa + (fic - 1)*this%dphi_fa
    end do

  end Subroutine
!----------------------------------------------------------------------------------
  function onestep_dtheodphi_nstep(this, rad, the, phi, ndiv) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the, phi
    integer, intent(in)::ndiv
    real*8::var
    integer fic
    real*8 dt
    var = the
    dt = phi/ndiv
    do fic = 1, ndiv
      var = onestep_dtheodphi(this, rad, var, dt)
    end do
  end function onestep_dtheodphi_nstep
!----------------------------------------------------------------------------------
  function onestep_dtheodphi(this, rad, the0, dt) result(var)
    ! dt is dphi
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, dt
    real*8::var
    real*8 the, dfdt, dfdt_sum

    the = the0
    !--1
    the = modulo(the, twopi)
    dfdt = 1.d0/this%getqloc_rt(rad, the)
    dfdt_sum = dfdt
    !--2
    the = the0 + dfdt*dt/2
    the = modulo(the, twopi)
    dfdt = 1.d0/this%getqloc_rt(rad, the)
    dfdt_sum = dfdt_sum + 2*dfdt
    !--3
    the = the0 + dfdt*dt/2
    the = modulo(the, twopi)
    dfdt = 1.d0/this%getqloc_rt(rad, the)
    dfdt_sum = dfdt_sum + 2*dfdt
    !--4
    the = the0 + dfdt*dt
    the = modulo(the, twopi)
    dfdt = 1.d0/this%getqloc_rt(rad, the)
    dfdt_sum = dfdt_sum + dfdt
    !--last
    var = the0 + dfdt_sum*dt/6
  end function onestep_dtheodphi

!----------------------------------------------------------------------------------
  subroutine equil_cls_readinput(this)
    implicit none
    class(equil_cls)::this
    integer iequmodel
    real*8 c1adhoc,c2adhoc,rmaxis_adhoc,Bmaxis_adhoc
    logical set_zerof  !--set fpol=0 manually
    character(100) fname
    character(100) profilename
    integer isource_ref
    real*8 rhoN, betaN, Bref, nref, Tref
    integer::nrad_sp, nthe_sp, nphi, nphiMult, nphifa_sp, nintg_fa
    real*8::radmax_sp, radmin_sp, nphieachfa

    integer r_order_bsp, z_order_bsp, phi_order_bsp, psi_order_bsp
    integer r_order_bsp3d, z_order_bsp3d, phi_order_bsp3d
    logical extrap
    integer isetpsi
    real*8 fBphi

    integer ifile, ierr, stat

!   integer,parameter::n_chars=1,n_real=1
!   character(100),dimension(n_chars)::char_parms
!   real*8::real_parms(n_real)

    namelist /equilibrium/ iequmodel,c1adhoc,c2adhoc,rmaxis_adhoc,Bmaxis_adhoc, &
      fname, profilename, isource_ref, Bref, nref, Tref, rhoN, betaN, set_zerof, &
      nrad_sp, nthe_sp, nphi, nphiMult, nphifa_sp, nintg_fa, &
      radmax_sp, radmin_sp, nphieachfa, &
      r_order_bsp, z_order_bsp, phi_order_bsp, psi_order_bsp, &
      r_order_bsp3d, z_order_bsp3d, phi_order_bsp3d, extrap, isetpsi, fBphi

    iequmodel=1
    c1adhoc=1.71
    c2adhoc=0.16
    rmaxis_adhoc=10.0
    Bmaxis_adhoc=2.d0

    set_zerof = .false.
    fname = 'g031213.00003'
    profilename = 'profile.h5'
    isource_ref = 1
    rhoN = 0.01
    betaN = 1e-6
    nref = 1d19
    Tref = 1000
    Bref = 1.d0

    r_order_bsp = 4
    z_order_bsp = 4
    phi_order_bsp = 4
    psi_order_bsp = 4

    r_order_bsp3d = 3
    z_order_bsp3d = 3
    phi_order_bsp3d = 3

    extrap = .true.

    isetpsi = 0
    fBphi = 1.d0

!   if (rank.eq.0) then
    ifile=201
    open (ifile, FILE='input',iostat=stat)
    if(stat.ne.0)write(*,*)' open input fails '
    read (ifile, equilibrium)
    close (ifile)
    if (rank .eq. 0) write (*, equilibrium)
!    char_parms(1)(1:100)=fname(1:100)
!    real_parms(1)=rhoN
!   endif

!   call mpi_bcast(char_parms,n_chars*100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!   call mpi_bcast(real_parms,n_real,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

!   if(rank/=0) then
!    fname(1:100)=char_parms(1)(1:100)
!    rhoN=real_parms(1)
!    write(*,equilibrium)
!   endif

    this%iequmodel=iequmodel
    !----adhoc----
    this%c1adhoc=c1adhoc
    this%c2adhoc=c2adhoc
    this%rmaxis_adhoc=rmaxis_adhoc
    this%Bmaxis_adhoc=Bmaxis_adhoc
    !----EQDSK----
    this%fname = fname
    this%profilename = profilename
    this%isource_ref = isource_ref
    this%Bref = Bref
    this%nref = nref
    this%Tref = Tref
    this%rhoN = rhoN
    this%betaN = betaN
    this%set_zerof = set_zerof
    !--set Spline&grid Parms
    this%nrad_sp = nrad_sp
    this%nthe_sp = nthe_sp
    this%nphi = nphi
    this%nphiMult = nphiMult
    this%nphifa_sp = nphifa_sp
    this%nintg_fa = nintg_fa

    this%radmax_sp = radmax_sp
    this%radmin_sp = radmin_sp
    this%nphieachfa = nphieachfa

    !--set internally
    this%extrap = .true.
!   this%extrap=.false.
    this%r_order_bsp = r_order_bsp
    this%z_order_bsp = z_order_bsp
    this%phi_order_bsp = phi_order_bsp
    this%psi_order_bsp = psi_order_bsp

    this%r_order_bsp3d = r_order_bsp3d
    this%z_order_bsp3d = z_order_bsp3d
    this%phi_order_bsp3d = phi_order_bsp3d

    this%extrap = extrap

    this%fBphi = fBphi
    this%isetpsi = isetpsi
  end Subroutine equil_cls_readinput

!----------------------------------------------------------------------------------
  subroutine equil_cls_calc_derived(this)
    implicit none
    class(equil_cls)::this
    real*8 vA, vth

    if (rank .eq. 0) write (*, *) '>>>>CALC_DERIVED,ISOURCE_REF=', this%isource_ref
    if (this%isource_ref .eq. 1) then
      this%Tref = (this%rhoN*charge_unit*this%Bref)**2/(2*mass_unit*charge_unit)
      this%nref = this%betaN*this%Bref**2/(2.d0*this%Tref*charge_unit*mu0)
      !--Extra Info
      vA = this%Bref/sqrt(mu0*this%nref*mass_unit)
      vth = sqrt(2.d0*this%Tref*charge_unit/mass_unit)
    elseif (this%isource_ref .eq. 2) then
      vA = this%Bref/sqrt(mu0*this%nref*mass_unit)
      vth = sqrt(2.d0*this%Tref*charge_unit/mass_unit)

      this%rhoN = sqrt(2.d0*mass_unit*this%Tref*charge_unit)/(charge_unit*this%Bref)
      this%betaN = vth**2/vA**2  !BETA=2T*c_u*mu0*n/B^2 ->n=beta*B^2/(2T*c_u*mu0)
    end if
    if(rank.eq.0)write(*,'(A6,e12.5,A6,e12.5,A6,e12.5,A6,e12.5,A6,e12.5,A6,e12.5)')'>>vA=',vA,',vth=',vth,',rhoN=',this%rhoN,  &
      ',betaN=', this%betaN, ',nref=', this%nref, ',Tref=', this%Tref
  end Subroutine equil_cls_calc_derived
!----------------------------------------------------------------------------------
  subroutine equil_cls_showinfo(this)
    implicit none
    class(equil_cls)::this

    character(LEN=80)::sform_i, sform_r

    if (rank .eq. 0) then
      sform_r = "(A14,e12.5,A10,e12.5,A10,e12.5,A10,e12.5,A10,e12.5,A10,e12.5)"
      sform_i = "(A14,I10,A14,I10,A14,I10,A14,I10,A14,I10,A14,I10)"
      write (*, sform_i) 'iequmodel=',this%iequmodel
      write (*, sform_r) 'c1adhoc=',this%c1adhoc,',c2adhoc=',this%c2adhoc,',rmaxis_adhoc=',this%rmaxis_adhoc,',Bmaxis_adhoc=',this%Bmaxis_adhoc
      write (*, sform_i) 'isource_ref=', this%isource_ref
      write (*, sform_r) 'rhoN=', this%rhoN, ',betaN=', this%betaN, ',nref=', this%nref, ',Tref=', this%Tref, ',Bref=', this%Bref
      write (*, sform_i) 'isetpsi=', this%isetpsi
      write (*, sform_r) 'fBphi=', this%fBphi
    end if
  end Subroutine equil_cls_showinfo
!----------------------------------------------------------------------------------
  function getradRZ(this, RR, ZZ) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::RR, ZZ
    real*8 var
    integer idx, idy, iflag

   if(this%iequmodel.eq.1)then
    idx = 0
    idy = 0
    call this%radrz_bsp%evaluate(RR, ZZ, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getradRZ:', iflag
   elseif(this%iequmodel.eq.2)then
    var=sqrt((RR-this%rmaxis)**2+(ZZ-this%zmaxis)**2)
   endif
  end Function
!----------------------------------------------------------------------------------
  function gettheRZ(this, RR, ZZ) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::RR, ZZ
    real*8 var
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    var = atan2(ZZ - this%zmaxis, RR - this%rmaxis)
    var = modulo(var, twopi)
   elseif(this%iequmodel.eq.2)then
    rormaj=this%getradRZ(RR,ZZ)/this%rmaxis
    var = atan2(ZZ - this%zmaxis, RR - this%rmaxis)
    var=2.d0*atan( sqrt((1.d0-rormaj)/(1.d0+rormaj))*tan(var/2.d0) )
    var = modulo(var, twopi)
    if(.not.(var.ge.0.d0.and.var.le.twopi))then
      var=0.d0
      write(*,*)'warning: var reset to 0 in gettheRZ ieqmodel==2 from ',var
    endif
   endif
  end Function
!----------------------------------------------------------------------------------
  subroutine calcBJgij(this, rad, the0, RR, ZZ, BB, dBdrad, dBdthe, Bthe_ct, Bphi_ct, jaco2, FF, g11, g12, g22)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8::RR, ZZ, BB, dBdrad, dBdthe, Bthe_ct, Bphi_ct, jaco2, FF, g11, g12, g22
    real*8 the, drdR, drdZ, dtdR, dtdZ

    the = modulo(the0, twopi)
    call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    g11 = drdR**2 + drdZ**2 !fix drdR->drdZ 230418
    g22 = dtdR**2 + dtdZ**2 !fix drdR->drdZ 230418
    g12 = drdR*dtdR + drdZ*dtdZ

    FF = this%getfpol_r(rad)
    RR = this%getR(rad, the)
    ZZ = this%getZ(rad, the)

    !BB=sqrt( (this%getdpsidr(rad)**2*g11+FF**2)/RR**2 )
    BB = this%getB(rad, the)
    dBdrad = this%getdBdrad(rad, the)
    dBdthe = this%getdBdthe(rad, the)

   if(this%iequmodel.eq.1)then
    Bthe_ct = this%getdpsidr(rad)/(jaco2*RR)
    Bphi_ct = FF/RR**2
   elseif(this%iequmodel.eq.2)then
    Bthe_ct=this%getBthe_ct(rad,the0)
    Bphi_ct = FF/RR**2
    !Bphi_ct=this%getBphi_ct(rad,the0)
   endif
  end Subroutine
!----------------------------------------------------------------------------------
  function geth_adhoc( rad, the, rmaxis) result(var)
    implicit none
    real*8 rad,the,rmaxis
    real*8 var
    real*8 rormaj

    rormaj=rad/rmaxis
    var=(1.d0-rormaj**2)/(1.d0-rormaj*cos(the))
  end Function
!----------------------------------------------------------------------------------
  function getBthe_ct(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    real*8 the, drdR, drdZ, dtdR, dtdZ, jaco2, RR

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    RR = this%getR(rad, the)
    var = this%getdpsidr(rad)/(jaco2*RR)
   elseif(this%iequmodel.eq.2)then
    var=this%Bmaxis/(this%rmaxis*this%getqloc_rt(rad,the0)*geth_adhoc(rad,the0,this%rmaxis)**2)
   endif
  end Function
!----------------------------------------------------------------------------------
  function getBphi_ct(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    var = this%getfpol_r(rad)/this%getR(rad, the)**2
   elseif(this%iequmodel.eq.2)then
    var=this%Bmaxis/(this%rmaxis*geth_adhoc(rad,the0,this%rmaxis)**2)
   endif
  end Function
!----------------------------------------------------------------------------------
  function getfeq(this, rad, the0, icase) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    integer, intent(in)::icase
    real*8 var
    real*8 the
    the = modulo(the0, twopi)

    if (abs(icase) == 1) then
      var = 1.d0
    elseif (abs(icase) == 2 .or. abs(icase) == 12 .or. abs(icase) .eq. 22) then !for parallel derivative
      var = (this%getBphi_ct(rad, the)*this%Bref/this%getB(rad, the))**2  !fix Bref 20250507
    elseif (abs(icase) == 14) then !for Poisson
      var = (this%Bref/this%getB(rad, the))**2
    elseif (abs(icase) == 15) then !for adiabatic electrons
      var = 2.0d0 !due to T normalization
    elseif (abs(icase) == 16 .or. abs(icase) == 116) then  !make sure consistent with solver_ext_cls_getfeq3d
      var = 1.d0
    elseif (abs(icase) == 17 .or. abs(icase) == 27) then !for parallel gradient^2 in (r,theta,phi)
      var = (this%getBthe_ct(rad, the)/this%getB(rad, the))**2
    elseif (abs(icase) == 18 .or. abs(icase) == 28) then !for parallel gradient^2 in (r,theta,phi)
      var = this%getBthe_ct(rad, the)*this%getBphi_ct(rad, the)/this%getB(rad, the)**2
    elseif (abs(icase) == 19 .or. abs(icase) == 29) then !for parallel gradient^2 in (r,theta,phi)
      var = (this%getBphi_ct(rad, the)/this%getB(rad, the))**2
    elseif (abs(icase) == 31) then !for Ohm parallel gradient: theta
      var = this%getBthe_ct(rad, the)/this%getB(rad, the)
    elseif (abs(icase) == 32) then !for Ohm parallel gradient: phi
      var = this%getBphi_ct(rad, the)/this%getB(rad, the)
    else
      write (*, *) '****error in EQU%GETFEQ:wrong ICAE****'
      call exit
    end if
    if (icase .lt. 0) then
      var = -var
    end if
  end Function
!----------------------------------------------------------------------------------
  function getgij_fa(this, rad, the0, phirow, phicol, ij) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phirow, phicol
    integer, intent(in)::ij
    real*8 var
    real*8 the, drdR, drdZ, dtdR, dtdZ, jaco2
    real*8 dedr_row, dedr_col, dedt_row, dedt_col, dedp_row, dedp_col

    the = modulo(the0, twopi)
    if (ij == 0) then
      var = 1.d0
    elseif (ij == 1) then
      call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
      var = drdR**2 + drdZ**2
    elseif (ij == 2) then
      call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
      var = (drdR**2 + drdZ**2)*this.getdetadrad(rad, the, phicol) &
            + (drdR*dtdR + drdZ*dtdZ)*this.getdetadthe(rad, the, phicol)
    elseif (ij == 4) then
      call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
      var = (drdR**2 + drdZ**2)*this.getdetadrad(rad, the, phirow) &
            + (drdR*dtdR + drdZ*dtdZ)*this.getdetadthe(rad, the, phirow)
    elseif (ij == 3 .or. ij == 7) then
      var = 0.d0
    elseif (ij == 5) then
      call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
      dedr_row = this%getdetadrad(rad, the, phirow)
      dedr_col = this%getdetadrad(rad, the, phicol)
      dedt_row = this%getdetadthe(rad, the, phirow)
      dedt_col = this%getdetadthe(rad, the, phicol)
      dedp_row = this%getdetadphi(rad, the, phirow)
      dedp_col = this%getdetadphi(rad, the, phicol)
      var = dedr_col*dedr_row*(drdR**2 + drdZ**2) &
            + dedt_col*dedt_row*(dtdR**2 + dtdZ**2) &
            + dedp_col*dedp_row/this%getR(rad, the)**2 &
            + (dedr_row*dedt_col + dedr_col*dedt_row)*(drdR*dtdR + drdZ*dtdZ)
    elseif (ij == 6) then
      var = this%getdetadphi(rad, the, phirow)/this%getR(rad, the)**2
    elseif (ij == 8) then
      var = this%getdetadphi(rad, the, phicol)/this%getR(rad, the)**2
    elseif (ij == 9) then
      var = 1/this%getR(rad, the)**2
    else
      var = 0; 
      if (rank .eq. 0) write (*, *) '****Error in getgij:wrong ij****'
    end if
  end Function getgij_fa

!----------------------------------------------------------------------------------
  function getgij_rtp(this, rad, the0, ij) result(var)
! gij in (r,theta) coordinates
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    integer, intent(in)::ij
    real*8 var
    real*8 the, drdR, drdZ, dtdR, dtdZ, jaco2
    real*8 dedr_row, dedr_col, dedt_row, dedt_col, dedp_row, dedp_col

    the = modulo(the0, twopi)
    if (ij == 0) then
      var = 1.d0
    elseif (ij == 1) then
      call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
      var = drdR**2 + drdZ**2
    elseif (ij == 2 .or. ij == 4) then
      call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
      var = drdR*dtdR + drdZ*dtdZ
    elseif (ij == 3 .or. ij == 6 .or. ij == 7 .or. ij == 8) then
      var = 0.d0
    elseif (ij == 5) then
      call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
      var = dtdR**2 + dtdZ**2
    elseif (ij == 9) then
      var = 1/this%getR(rad, the)**2
    else
      var = 0; 
      if (rank .eq. 0) write (*, *) '****Error in getgij:wrong ij****'
    end if
  end Function getgij_rtp
!----------------------------------------------------------------------------------
  function getB(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 0
    call this%Brt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getB:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    var=this%Bmaxis*(1-rormaj*cos(the0))/(1-rormaj**2)*sqrt((rormaj/this%getqloc_rt(rad,the0))**2/(1.d0-rormaj**2)+1.d0)
   endif
  end Function
!----------------------------------------------------------------------------------
  function getdBdrad(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj,qloc
    real*8 tmpcoeff

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 1
    idy = 0
    call this%Brt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getB:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    qloc=this%getqloc_rt(rad,the0)
    tmpcoeff = 1-rormaj**2
    var=( ( 1.d0/rad-this%getdqdrloc_rt(rad,the0)/qloc  &
            +3*rormaj/(tmpcoeff)/this%rmaxis )*(rormaj/qloc)**2 &
            +2.d0*rormaj/this%rmaxis )*(1.d0-rormaj*cos(the0))*this%Bmaxis/(tmpcoeff)**2 &
       -( (rormaj/qloc)**2/(tmpcoeff)+1.d0 )*cos(the0)*this%Bmaxis/((tmpcoeff)*this%rmaxis)
    var=var/sqrt( rormaj**2/(qloc**2*(tmpcoeff))+1.d0 )
   endif
  end Function
!----------------------------------------------------------------------------------
  function getdBdthe(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj,qloc

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 1
    call this%Brt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getB:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    qloc=this%getqloc_rt(rad,the0)
    var=this%Bmaxis*rormaj*sin(the0)/(1.d0-rormaj**2)*sqrt( rormaj**2/(qloc**2*(1.d0-rormaj**2))+1.d0 )
   endif
  end Function
!----------------------------------------------------------------------------------
  function getBrad_co(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 0
    call this%Brad_co_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getBrad_co:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    !var=rormaj**2*this%Bmaxis*sin(the0)*(1.d0-rormaj*cos(the0))**2/( this%getqloc_rt(rad,the0)*sqrt(1-rormaj**2)**5 )
    var=rormaj**2*this%Bmaxis*sin(the0)/(this%getqloc_rt(rad,the0)*(1-rormaj**2)**2 )
   endif
  end Function
!----------------------------------------------------------------------------------
  function getdBrad_codrad(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 1
    idy = 0
    call this%Brad_co_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getdBrad_co/drad:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
!    var=( -rad*this%getdqdrloc_rt(rad,the0)/this%getqloc_rt(rad,the0)+2.d0 &
!          -2.d0*cos(the0)/(1-rormaj*cos(the0))*rormaj &
!          +5.d0*rormaj**2/(1.d0-rormaj**2) ) &
!        *rormaj*this%Bmaxis*sin(the0)*(1.d0-rormaj*cos(the0))**2/( this%getqloc_rt(rad,the0)*sqrt(1-rormaj**2)**5 ) !from Brad_co
    var=( -rad*this%getdqdrloc_rt(rad,the0)/this%getqloc_rt(rad,the0)+2.d0 &
          +4.d0*rormaj**2/(1.d0-rormaj**2) ) &
        *rormaj*this%Bmaxis*sin(the0)/( this%rmaxis*this%getqloc_rt(rad,the0)*(1-rormaj**2)**2 ) !from Brad_co
   endif
  end Function
!----------------------------------------------------------------------------------
  function getdBrad_codthe(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then

    the = modulo(the0, twopi)

    idx = 0
    idy = 1
    call this%Brad_co_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getdBrad_co/dthe:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    !var=( cos(the0)+2.d0*rormaj*sin(the0)**2/(1-rormaj*cos(the0)) ) &
    !    *rormaj**2*this%Bmaxis*(1.d0-rormaj*cos(the0))**2/( this%getqloc_rt(rad,the0)*sqrt(1-rormaj**2)**5 ) !from Brad_co
    var= cos(the0)*rormaj**2*this%Bmaxis/( this%getqloc_rt(rad,the0)*(1-rormaj**2)**2 )
        !*rormaj*this%Bmaxis*sin(the0)/( this%rmaxis**2*this%getqloc_rt(rad,the0)*sqrt(1-rormaj**2)**3 ) !from Brad_co
   endif
  end Function
!----------------------------------------------------------------------------------
  function getBthe_co(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 0
    call this%Bthe_co_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getBthe_co:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    !var=rormaj**2*this%Bmaxis*this%rmaxis*(1.d0-rormaj*cos(the0))**2/( this%getqloc_rt(rad,the0)*(1.d0-rormaj**2)**2 )
    var=rad**2*this%Bmaxis/( this%rmaxis*this%getqloc_rt(rad,the0)*(1.d0-rormaj**2) )
   endif
  end Function
!----------------------------------------------------------------------------------
  function getBphi_co(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    var = this%getfpol_r(rad)
   elseif(this%iequmodel.eq.2)then
    var = this%Bmaxis*this%rmaxis !/this%getR(rad,the0)
   endif
  end Function
!----------------------------------------------------------------------------------
  function getdBthe_codrad(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 1
    idy = 0
    call this%Bthe_co_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getdBthe_co/drad:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    !var=2.d0*( -rad*this%getdqdrloc_rt(rad,the0)/this%getqloc_rt(rad,the0)+1.d0 &
    !      -2.d0*cos(the0)/(1-rormaj*cos(the0))*rormaj &
    !      +2.d0*rormaj**2/(1.d0-rormaj**2) ) &
    !    *rormaj*this%Bmaxis*this%rmaxis*(1.d0-rormaj*cos(the0))**2/( this%getqloc_rt(rad,the0)*(1.d0-rormaj**2)**2 ) !from Bthe_co
    var=( -rad*this%getdqdrloc_rt(rad,the0)/this%getqloc_rt(rad,the0)+2.d0 &
          +2.d0*rormaj**2/(1.d0-rormaj**2) ) &
        *rormaj*this%Bmaxis/( this%getqloc_rt(rad,the0)*(1.d0-rormaj**2) )!from Bthe_co
   endif
  end Function
!----------------------------------------------------------------------------------
  function getdBthe_codthe(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 1
    call this%Bthe_co_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getdBthe_co/dthe:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    !var=2.d0*rormaj**3*this%Bmaxis*this%rmaxis*(1.d0-rormaj*cos(the0))*sin(the0) &
    !        /(this%getqloc_rt(rad,the0)*(1.d0-rormaj**2)**2)
    var=0.d0
   endif
  end Function

!----------------------------------------------------------------------------------
  function getBdirect(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    real*8 the, drdR, drdZ, dtdR, dtdZ, jaco2, FF, g11, RR

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    g11 = drdR**2 + drdZ**2 !fix drdR->drdZ 20230418

    FF = this%getfpol_r(rad)
    RR = this%getR(rad, the)
    var = sqrt((this%getdpsidr(rad)**2*g11 + FF**2)/RR**2)
   elseif(this%iequmodel.eq.2)then
    var=this%getB(rad,the0)
   endif
  end Function

!----------------------------------------------------------------------------------
  function getBrad_co_direct(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    real*8 the, drdR, drdZ, dtdR, dtdZ, jaco2, g12, RR

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    g12 = drdR*dtdR + drdZ*dtdZ

    RR = this%getR(rad, the)
    var = -jaco2*g12*this%getdpsidr(rad)/RR
   elseif(this%iequmodel.eq.2)then
    var=this%getBrad_co(rad,the0)
   endif
  end Function
!----------------------------------------------------------------------------------
  function getBthe_co_direct(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    real*8 the, drdR, drdZ, dtdR, dtdZ, jaco2, g11, RR

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    g11 = drdR**2 + drdZ**2

    RR = this%getR(rad, the)
    var = jaco2*g11*this%getdpsidr(rad)/RR
   elseif(this%iequmodel.eq.2)then
    var=this%getBthe_co(rad,the0)
   endif
  end Function
!----------------------------------------------------------------------------------
  subroutine calc_curlb_ct_direct(this, rad, the0, frad, fthe, fphi)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8, intent(out)::frad, fthe, fphi

    real*8 the, drdR, drdZ, dtdR, dtdZ, jaco2, jaco3, g11, g12, g22, RR
    real*8 Bphi_co, dBdrad, dBdthe, JB2inv, BB

    the = modulo(the0, twopi)

    call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    g11 = drdR**2 + drdZ**2
    g22 = dtdR**2 + dtdZ**2
    g12 = drdR*dtdR + drdZ*dtdZ

    RR = this%getR(rad, the)
    jaco3 = jaco2*RR
    Bphi_co = this%getfpol_r(rad)

    dBdrad = this%getdBdrad(rad, the)
    dBdthe = this%getdBdthe(rad, the)
    BB = this%getB(rad, the)
    JB2inv = 1.d0/(jaco3*BB**2)

    frad = Bphi_co*dBdthe*JB2inv
    fthe = -Bphi_co*dBdrad*JB2inv + (g11*g22 - g12**2)*jaco3*this%getdfpoldr_r(rad)/(BB*RR**2)
    fphi = -(this%getBrad_co(rad, the)*dBdthe - this%getBthe_co(rad, the)*dBdrad)*JB2inv &
           + (this%getdBrad_codthe(rad, the) - this%getdBthe_codrad(rad, the))/(jaco3*BB)
  end Subroutine Calc_curlb_ct_direct
!----------------------------------------------------------------------------------
  function getcurlbrad_ct(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 0
    call this%curlbrad_ct_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getcurlbrad_ct:', iflag
   elseif(this%iequmodel.eq.2)then
    var=this%getBphi_co(rad,the0)*this%getdBdthe(rad,the0)/(this%getjaco3(rad,the0)*this%getB(rad,the0)**2)
   endif
  end Function
!----------------------------------------------------------------------------------
  function getcurlbthe_ct(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 0
    call this%curlbthe_ct_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getcurlbthe_ct:', iflag
   elseif(this%iequmodel.eq.2)then
    var= -this%getBphi_co(rad,the0)*this%getdBdrad(rad,the0)/(this%getjaco3(rad,the0)*this%getB(rad,the0)**2)
   endif
  end Function
!----------------------------------------------------------------------------------
  function getcurlbphi_ct(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 jaco3,BB

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 0
    call this%curlbphi_ct_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (rank .eq. 0 .and. iflag .ne. 0) write (*, *) '****Error in getcurlbphi_ct:', iflag
   elseif(this%iequmodel.eq.2)then
    BB=this%getB(rad,the0)
    jaco3=this%getjaco3(rad,the0)
    var= -( this%getBrad_co(rad, the)*this%getdBdthe(rad,the0) - this%getBthe_co(rad, the)*this%getdBdrad(rad,the0) ) &
          /(jaco3*BB**2) &
         + (this%getdBrad_codthe(rad, the) - this%getdBthe_codrad(rad, the))/(jaco3*BB)
   endif
  end Function

!----------------------------------------------------------------------------------
  function getthe3d(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 0
    idy = 0
    idz = 0
    call this%the3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in getthe3d:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=the0+phi/this%getqloc_rt(rad,the0)
   endif
  end
!----------------------------------------------------------------------------------
  function getdthe3ddrad(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 1
    idy = 0
    idz = 0
    call this%the3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in getdthe3ddrad:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=-phi*this%getdqdrloc_rt(rad,the0)/this%getqloc_rt(rad,the0)**2
   endif
  end
!----------------------------------------------------------------------------------
  function getdthe3ddthe(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 0
    idy = 1
    idz = 0
    call this%the3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in getdthe3ddthe:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=1.d0
   endif
  end
!----------------------------------------------------------------------------------
  function getdthe3ddphi(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 0
    idy = 0
    idz = 1
    call this%the3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in getdthe3ddphi:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=1.d0/this%getqloc_rt(rad,the0)
   endif
  end
!----------------------------------------------------------------------------------
  function geteta(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 0
    idy = 0
    idz = 0
    call this%eta3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in geteta:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=the0-phi/this%getqloc_rt(rad,the0)
   endif

  end
!----------------------------------------------------------------------------------
  function getdetadrad(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 1
    idy = 0
    idz = 0
    call this%eta3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in getdetadrad:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=+phi*this%getdqdrloc_rt(rad,the0)/this%getqloc_rt(rad,the0)**2
   endif
  end
!----------------------------------------------------------------------------------
  function getdetadthe(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 0
    idy = 1
    idz = 0
    call this%eta3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in getdetadthe:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=1.d0
   endif
  end
!----------------------------------------------------------------------------------
  function getdetadphi(this, rad, the0, phi) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0, phi
    real*8 var
    integer idx, idy, idz, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 0
    idy = 0
    idz = 1
    call this%eta3d_bsp%evaluate(rad, the, phi, idx, idy, idz, var, iflag)
    if (iflag .ne. 0) write (*, '(A30,I6,3e14.5)') '****Error in getdetadphi:', iflag, rad, the, phi
   elseif(this%iequmodel.eq.2)then
    var=-1.d0/this%getqloc_rt(rad,the0)
   endif
  end
!----------------------------------------------------------------------------------
  function getcos_straight_adhoc(rad, the0,rmaxis) result(var)
    implicit none
    real*8, intent(in)::rad, the0,rmaxis
    real*8 var
    real*8 rormaj

    rormaj=rad/rmaxis
    var= -1.d0/rormaj+(1-rormaj**2)/(rormaj*(1.d0-rormaj*cos(the0)))
    !var= (cos(the0)-rormaj)/(1.d0-rormaj*cos(the0))
  end
!----------------------------------------------------------------------------------
  function getsin_straight_adhoc( rad, the0,rmaxis) result(var)
    implicit none
    real*8, intent(in)::rad, the0,rmaxis
    real*8 var
    real*8 rormaj

    rormaj=rad/rmaxis
    var= sqrt(1-rormaj**2)*sin(the0)/(1.d0-rormaj*cos(the0))
  end
!----------------------------------------------------------------------------------
  function getR(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)
    idx = 0
    idy = 0
    call this%rrt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (iflag .ne. 0) write (*, *) '****Error in getR:', iflag
   elseif(this%iequmodel.eq.2)then
    !var=this%rmaxis+rad*cos(the0)
    var=this%rmaxis+rad*getcos_straight_adhoc(rad,the0,this%rmaxis)
   endif
  end

!----------------------------------------------------------------------------------
  function getdRdrad(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

    the = modulo(the0, twopi)
   if(this%iequmodel.eq.1)then

    idx = 1
    idy = 0
    call this%rrt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (iflag .ne. 0) write (*, *) '****Error in getdRdrad:', iflag
   elseif(this%iequmodel.eq.2)then
    !var=cos(the0)
    rormaj=rad/this%rmaxis
    var=getcos_straight_adhoc(rad,the0,this%rmaxis) &
            -rormaj*sin(the0)**2/(1.d0-rormaj*cos(the0))**2
   endif
  end
!----------------------------------------------------------------------------------
  function getdRdthe(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 1
    call this%rrt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (iflag .ne. 0) write (*, *) '****Error in getdRdthe:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    !var=-rad*sin(the0)
    var=-rad*(1.d0-rormaj**2)*sin(the0)/(1.d0-rormaj*cos(the0))**2
   endif
  end
!----------------------------------------------------------------------------------
  function getZ(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 0
    call this%zrt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (iflag .ne. 0) write (*, *) '****Error in getZ:', iflag
   elseif(this%iequmodel.eq.2)then
    !var=this%zmaxis+rad*sin(the0)
    var=this%zmaxis+rad*getsin_straight_adhoc(rad,the0,this%rmaxis)
   endif
  end
!----------------------------------------------------------------------------------
  function getdZdrad(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 1
    idy = 0
    call this%zrt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (iflag .ne. 0) write (*, *) '****Error in getdZdrad:', iflag
   elseif(this%iequmodel.eq.2)then
    !var=sin(the0)
    rormaj=rad/this%rmaxis
    var=getsin_straight_adhoc(rad,the0,this%rmaxis) &
            +rormaj*sin(the0)*(cos(the0)-rormaj)/(sqrt(1.d0-rormaj**2)*(1-rormaj*cos(the0))**2)
   endif
  end
!----------------------------------------------------------------------------------
  function getdZdthe(this, rad, the0) result(var)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the0
    real*8 var
    integer idx, idy, iflag
    real*8 the
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    the = modulo(the0, twopi)

    idx = 0
    idy = 1
    call this%zrt_bsp%evaluate(rad, the, idx, idy, var, iflag)
    if (iflag .ne. 0) write (*, *) '****Error in getdZdthe:', iflag
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    !var=rad*cos(the0)
    var=rad*sqrt(1.d0-rormaj**2)/(1.d0-rormaj*cos(the0))**2*(cos(the0)-rormaj)
   endif
  end
!----------------------------------------------------------------------------------
  function getdqdrloc_rt(this, rad, the) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad, the
    real*8 var
    real*8 drdR, drdZ, dtdR, dtdZ, jaco2, RR
    real*8 rormaj

   if(this%iequmodel.eq.1)then
    write(*,*)'error: not implemented'
    call exit
   elseif(this%iequmodel.eq.2)then
    rormaj=rad/this%rmaxis
    var= 2.d0*this%c2adhoc*rad &
        +(this%c1adhoc+this%c2adhoc*rad**2)*rormaj/( (1.d0-rormaj**2)*this%rmaxis )
    var = var /sqrt(1.d0-rormaj**2)
   endif
  end
!----------------------------------------------------------------------------------
  function getqloc_rt(this, rad, the) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad, the
    real*8 var
    real*8 drdR, drdZ, dtdR, dtdZ, jaco2, RR

   if(this%iequmodel.eq.1)then
    call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    RR = this%getR(rad, the)
    var = jaco2*this%getfpol_r(rad)/(RR*this%getdpsidr(rad))
   elseif(this%iequmodel.eq.2)then
    var=this%c1adhoc+this%c2adhoc*rad**2
    var=var/sqrt(1.d0-(rad/this%rmaxis)**2)
   endif
  end
!----------------------------------------------------------------------------------
  function getpsi(this, rad) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad
    real*8 var

   if(this%iequmodel.eq.1)then
    var = this%simag + this%psiwid*rad**2
   elseif(this%iequmodel.eq.2)then
    var=this%Bmaxis/(2*this%c2adhoc)*log(1.d0+this%c2adhoc/this%c1adhoc*rad**2)
   endif
  end
!----------------------------------------------------------------------------------
  function getdpsidr(this, rad) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad
    real*8 var

   if(this%iequmodel.eq.1)then
    var = 2*this%psiwid*rad
   elseif(this%iequmodel.eq.2)then
    var=this%Bmaxis*rad/(this%c1adhoc+this%c2adhoc*rad**2)
   endif
  end
!----------------------------------------------------------------------------------
  function getfpol_r(this, rad) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad
    real*8 var
    integer idx, iflag

   if(this%iequmodel.eq.1)then
    idx = 0
    call this%F1d_bsp%evaluate(this%getpsi(rad), idx, var, iflag)
   if (iflag .ne. 0) write (*, '(A20,I5,A5,e15.5,A5,e15.5)') '****Error in getfpol_r:', iflag, ',r=', rad, ',psi=', this%getpsi(rad)
   elseif(this%iequmodel.eq.2)then
    var=this%Bmaxis*this%rmaxis
   endif
  end
!----------------------------------------------------------------------------------
  function getdfpoldr_r(this, rad) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad
    real*8 var
    integer idx, iflag

   if(this%iequmodel.eq.1)then
    idx = 1
    call this%F1d_bsp%evaluate(this%getpsi(rad), idx, var, iflag)
    var = var*2*this%psiwid*rad
if (iflag .ne. 0) write (*, '(A20,I5,A5,e15.5,A5,e15.5)') '****Error in getdfpoldr_r:', iflag, ',r=', rad, ',psi=', this%getpsi(rad)
   elseif(this%iequmodel.eq.2)then
    var=0.d0
   endif
  end
!----------------------------------------------------------------------------------
  function getjaco2(this, rad, the) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad, the
    real*8 var
    real*8 dRdrad, dRdthe, dZdrad, dZdthe
    integer idx, idy, iflag

    dRdrad = this%getdRdrad(rad, the)
    dRdthe = this%getdRdthe(rad, the)
    dZdrad = this%getdZdrad(rad, the)
    dZdthe = this%getdZdthe(rad, the)

    var = dRdrad*dZdthe - dRdthe*dZdrad
  end
!----------------------------------------------------------------------------------
  function getjaco3(this, rad, the) result(var)
    implicit none
    class(equil_cls)::this
    real*8 rad, the
    real*8 var
    real*8 dRdrad, dRdthe, dZdrad, dZdthe, RR
    integer idx, idy, iflag

    dRdrad = this%getdRdrad(rad, the)
    dRdthe = this%getdRdthe(rad, the)
    dZdrad = this%getdZdrad(rad, the)
    dZdthe = this%getdZdthe(rad, the)
    RR = this%getR(rad, the)
    var = RR*(dRdrad*dZdthe - dRdthe*dZdrad)
  end
!----------------------------------------------------------------------------------
  subroutine calcdrtdRZ(this, rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
    implicit none
    class(equil_cls)::this
    real*8, intent(in)::rad, the
    real*8, intent(out)::drdR, drdZ, dtdR, dtdZ, jaco2
    real*8 dRdrad, dRdthe, dZdrad, dZdthe, jacoinv
    integer idx, idy, iflag

    dRdrad = this%getdRdrad(rad, the)
    dRdthe = this%getdRdthe(rad, the)
    dZdrad = this%getdZdrad(rad, the)
    dZdthe = this%getdZdthe(rad, the)

    jaco2 = dRdrad*dZdthe - dRdthe*dZdrad
    !jaco2 = rad*sqrt(1-rad**2/this%rmaxis**2)/(1-rad/this%rmaxis*cos(the)) !dRdrad*dZdthe - dRdthe*dZdrad
    jacoinv = 1/jaco2
    drdR = dZdthe*jacoinv
    drdZ = -dRdthe*jacoinv
    dtdR = -dZdrad*jacoinv
    dtdZ = dRdrad*jacoinv
  end
!----------------------------------------------------------------------------------
  subroutine record_var2d(this)
    implicit none
    class(equil_cls)::this
    integer fic, fjc, nfile,stat
    real*8 rad, the, RR, ZZ, drdR, drdZ, dtdR, dtdZ, jaco2, BB, dBdrad, dBdthe, Bthe_ct, Bphi_ct, FF, g11, g12, g22
    real*8 curlbrad_ct, curlbthe_ct, curlbphi_ct, Brad_co, Bthe_co, Bphi_co, Babs

    if (rank .eq. 0) then
      write(*,*)'----Equilibrium record_var2d starts----'
      nfile = 201
      open (nfile, file='equvar2d.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open equvar2d fails '
      do fjc = 1, this%nthe_sp
        do fic = 1, this%nrad_sp
          rad = this%rad1d_sp(fic)
          the = this%the1d_sp(fjc)
          call this%calcdrtdRZ(rad, the, drdR, drdZ, dtdR, dtdZ, jaco2)
          call this%calcBJgij(rad, the, RR, ZZ, BB, dBdrad, dBdthe, Bthe_ct, Bphi_ct, jaco2, FF, g11, g12, g22)
          !call this%calc_curlb_ct_direct(rad,the,curlbrad_ct,curlbthe_ct,curlbphi_ct)
          curlbrad_ct = this%getcurlbrad_ct(rad, the)
          curlbthe_ct = this%getcurlbthe_ct(rad, the)
          curlbphi_ct = this%getcurlbphi_ct(rad, the)
          Brad_co = this%getBrad_co(rad, the)
          Bthe_co = this%getBthe_co(rad, the)
          Bphi_co = this%getBphi_co(rad, the)
          Babs = sqrt(Bthe_co*Bthe_ct + Bphi_co*Bphi_ct)
          !1-5:rad,the,R,Z,q, 6-10:drdR,drdZ,dtdR,dtdZ,j2,11:B,dBdrad,dBdthe,Bthe_ct,Bphi_ct,
          !16:fpol,g11,g12,g22,curbrad_ct,21:curbthe_ct,curlbphi_ct,Brad_co,Bthe_co,Babs
          write (nfile, '(50e16.4)') rad, the, RR, ZZ, this%getqloc_rt(rad, the), &
            drdR, drdZ, dtdR, dtdZ, jaco2, &
            BB, dBdrad, dBdthe, Bthe_ct, Bphi_ct, &
            FF, g11, g12, g22, curlbrad_ct, &
            curlbthe_ct, curlbphi_ct, Brad_co, Bthe_co, Bphi_co, &
            Babs
        end do
      end do
      close (nfile)
      write(*,*)'----Equilibrium record_var2d ends----'
    end if !RANK==0
  end
!----------------------------------------------------------------------------------
  subroutine record_var3d(this)
    implicit none
    class(equil_cls)::this
    integer fic, fjc, fkc, nfile, stat
    real*8 rad, the, phi, eta, RR, ZZ

    if (rank .eq. 0) then
      write(*,*)'----Equilibrium record_var3d starts----'
      nfile = 201
      open (nfile, file='eta3d.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open eta3d.txt fails '
      do fkc = 1, this%nphifa_sp
        do fjc = 1, this%nthe_sp
          do fic = 1, this%nrad_sp
            rad = this%rad1d_sp(fic)
            the = this%the1d_sp(fjc)
            phi = this%phifa1d_sp(fkc)
            eta = this%geteta(rad, the, phi)
            RR = this%getR(rad, eta)
            ZZ = this%getZ(rad, eta)
            !R,Z,eta,deta/drad,deta/dthe,deta/dphi
         write(nfile,'(20e16.4)')RR,ZZ,eta,this%getdetadrad(rad,the,phi),this%getdetadthe(rad,the,phi),this%getdetadphi(rad,the,phi)
          end do
        end do
      end do
      close (nfile)
    end if !RANK==0

    if (rank .eq. mpisize - 1) then
      nfile = 202
      open (nfile, file='the3d.txt',iostat=stat)
      if(stat.ne.0)write(*,*)' open the3d.txt fails '
      do fkc = 1, this%nphifa_sp
        do fjc = 1, this%nthe_sp
          do fic = 1, this%nrad_sp
            ! eta is used as theta;the as eta
            rad = this%rad1d_sp(fic)
            eta = this%the1d_sp(fjc)
            phi = this%phifa1d_sp(fkc)
            the = this%getthe3d(rad, eta, phi)
            RR = this%getR(rad, the)
            ZZ = this%getZ(rad, the)
            !R,Z,eta,dthe/drad,dthe/deta,dthe/dphi
   write(nfile,'(20e16.4)')RR,ZZ,the,this%getdthe3ddrad(rad,eta,phi),this%getdthe3ddthe(rad,eta,phi),this%getdthe3ddphi(rad,eta,phi)
          end do
        end do
      end do
      close (nfile)
      write(*,*)'----Equilibrium record_var3d ----'
    end if !RANK==0
  end
!----------------------------------------------------------------------------------
end module
!----------------------------------------------------------------------------------
