c=======================================================================
c GCOM-C OceanColor atmospheric correction
c This code is used with ../bin/run_pub_sgli_Global.sh, ../lut
c see ../readme.txt
c
c Copyright (c) Hiroshi Murakami, Japan Aerospace Exploration Agency (JAXA)
c 
c Permission is hereby granted, free of charge, to any person obtaining a 
c copy of this software and associated documentation files (the "Software"), 
c to deal in the Software without restriction, including without limitation
c the rights to use, copy, modify, merge, publish, distribute, sublicense, 
c and/or sell copies of the Software, subject to the following conditions:
c 1. The Software is provided "as is", without warranty of any kind, express
c   or implied, including but not limited to the warranties of merchantability,
c   fitness for a particular purpose and noninfringement. In no event shall the
c   authors or copyright holders be liable for any claim, damages or other 
c   liability, whether in an action of contract, tort or otherwise, arising 
c   from, out of or in connection with the Software or the use or other dealings
c   in the Software.
c 2. The authors reserve the right to change the terms of use at any time without notice.
c 3. The authors reserve the right to discontinue distribution of the Software 
c   at any time without notice.
c 4. The use of the Software for any purpose that is illegal, immoral, or 
c   otherwise deemed inappropriate by societal standards is prohibited.
c By using the Software, you are consenting to the above terms.
c=======================================================================
c
      Program pub_sgli_Global
c
      USE HDF5	! This module contains all necessary modules
c
      implicit none
c
      external read_tbl, interp_tbl,read_ANC, read_TOZ, SUNGLINT, 
     & estimate_aot, rw_sim_lci, brdf_oc, cal_iop, cal_par, read2d_full,
     & h5read_r8, read_dattr, h5_write, h5_write_r8, gl_attr, sp_attr
c
c     ------------------------------------------------------------------
      character(LEN=4)	:: avsn=  '2.09' ! algorithm version (0~z)
      character(LEN=6)	:: pvsn='009.11' ! parameter version (0000~zzzz)
      character(LEN=4)	:: fvsn	! product version (avsn(1:1)//pvsn(1:3))
      integer,parameter :: nmdl=9		! 10:v1, 9*:v2
      integer			:: ifg_brdf=1	! 0: off, 1*: on
      integer			:: ifg_sung=1	! 0: off, 1*: on
      integer			:: ifg_b2ss=1	! 0: off, 1*: on (bbp->tsm, adg->ag)
      integer			:: ifngc=1		! 1*:neg rrs corr on, 0:off
      integer			:: ifg_rvsc=1	! 0: off, 1*: on
      integer			:: ifluce=1		! 0: off, 1*:LUCA error corr on
c     ----------------------------------------------------
c
      real(4),parameter :: pi=3.14159265,d2r=pi/180.
      real(4),parameter	:: tsg1=0.16	! sunglint mask (see attr of QA)
      real(4),parameter	:: tsg2=0.01	! sunglint flag (see attr of QA)
      real(4),parameter	:: rcmin=0.0001	! min reflectance
      integer,parameter :: iw0=2	! window for cloud flag
c
      integer,parameter :: nbl0=15
      real(4) :: wvl(nbl0)
     &       =(/ 380.0296, 412.5107, 443.2402, 489.8488, 529.6403,
     &           566.1529, 671.9961, 672.0980, 763.0741, 866.7649,
     &           867.1204,1054.9942,1385.3511,1634.5061,2209.4816/)
      real(4) :: f0n(nbl0)
     &       =(/1092.1439,1712.1658,1898.3187,1938.4599,1850.9603,
     &          1797.1381,1502.5498,1502.3006,1245.4534, 956.3376,
     &           956.6228, 646.5423, 361.2380, 237.5831,  84.2457/)
      real	:: bw(nbl0)=
     &      (/ 10.7, 10.4, 10.2, 10.4, 19.7,
     &         20.1, 22.3, 22.1, 11.4, 21.1,
     &         21.3, 21.4, 20.2,195.8, 50.8/)
c
      real(4) :: kv(nbl0),kt(nbl0),kb(nbl0)
      real(4) :: kv0(nbl0)=(/
     & 0.987,1.034,1.014,1.027,1.060,1.040,1.000,0.991,1.000,1.000,
     & 1.000,1.000,1.000,1.000,1.000/)	! MOBY*BOUSSOLE (adj Day 950)
c
      real(4) :: kt0(nbl0)=(/
     &-6.20E-5,-6.06E-5,-5.79E-5,-5.05E-5,-4.20E-5,-3.06E-5,-5.04E-6,
     &-4.34E-6, 0., 0., 0., -4.20E-5,-5.82E-5, 0., 0./)
c
      real(4) :: kbq(3,nbl0)=reshape( (/
     & -0.021426294,0.001064101,-2.48E-07,
     &  0.008461326,0.000747790,-4.41E-08,
     & -0.002959722,0.000607677, 8.68E-10,
     & -0.010513429,0.000653283,-1.21E-07,
     &  0.003113079,0.000237091, 1.66E-08,
     & -0.012615250,0.000408857,-1.20E-07,
     & -0.011459627,0.000252642,-4.15E-08,
     & -0.004370468,0.000187907,-1.21E-08,
     &  0.008436488,0.000222141, 3.16E-08,
     & -0.000424147,0.000274879,-1.17E-07,
     & -0.000248941,0.000219977,-1.67E-08,
     &  0.,0.,0., 0.,0.,0., 0.,0.,0., 0.,0.,0.
     &  /),(/3,nbl0/) )
      real(4) :: kbk(3,nbl0)=reshape( (/
     & -0.009722052,0.000835637,-1.81E-07,
     &  0.007593721,0.000688787,-1.71E-07,
     & -0.001504649,0.000536011,-6.42E-08,
     & -0.002564131,0.000509960,-5.46E-08,
     &  0.003147961,0.000215445,-2.95E-08,
     & -0.003940712,0.000303020,-5.33E-08,
     & -0.005830028,0.000197740,-1.84E-08,
     & -0.003119111,0.000169235,-4.63E-08,
     &  0.006613161,0.000194620,-3.59E-08,
     & -0.003417236,0.000269736,-1.52E-07,
     & -0.000727193,0.000195501,-5.33E-08,
     &  0.,0.,0., 0.,0.,0., 0.,0.,0., 0.,0.,0.
     &  /),(/3,nbl0/) )
c
      real(4) :: kte(11,4)=reshape( (/
     & 9.984e-01, 1.001e+00, 1.000e+00, 1.001e+00, 1.001e+00, 1.002e+00,
     & 1.001e+00, 1.003e+00, 1.000e+00, 9.994e-01, 9.989e-01,
     & 1.190e-05,-2.898e-06,-1.840e-06,-5.527e-06,-7.283e-06,-7.996e-06,
     &-3.405e-06,-1.459e-05, 2.368e-07, 3.307e-06, 4.889e-06,
     &-1.881e-08,-7.158e-10, 9.951e-10, 5.597e-09, 6.358e-09, 6.333e-09,
     &-1.124e-09, 1.178e-08,-5.283e-09,-7.121e-09,-6.329e-09,
     & 6.443e-12, 1.113e-12, 7.615e-13,-9.625e-13,-9.631e-13,-7.691e-13,
     & 1.189e-12,-2.548e-12, 2.038e-12, 2.260e-12, 1.535e-12
     & /),(/11,4/) )
c
      real(4) :: lcre(11)
c
      REAL(4) :: rpx1(201),rpx2(201)	! RVS (n-1.)/(nl-1.)*200.+1
      REAL(4) :: rp,rpix(5000,nbl0)
c
      integer,parameter :: nbul=11+2
      integer :: nbu(nbul)=(/1,2,3,4,5,6,7,8,9,10,11,12,14/)
      real(4) :: slpr(nbul),offr(nbul)
      real(4) :: dra,dse
      real(4) :: f0(nbul)
      real(4) :: wl(nbul)
      integer,parameter :: ngl=7
      real(4) :: slpg(ngl),offg(ngl)
      character(LEN=4)	:: l1vsn
      integer,parameter :: nbal=2
      integer :: nba(nbal)=(/10,7/)	! in nbu
c
      integer,parameter :: nbwl0=7	! VN01~07
      character(LEN=3)	:: nbw(nbwl0)=(/'380','412','443','490','530',
     &                                  '565','670'/)
      real(4) :: slpw(nbwl0),offw(nbwl0)
c
c     PAR
      integer,parameter :: nbpl=6
      integer :: nbpi(nbpl)=(/2,3,5,8,11,12/)	! in nbu
      real(4) :: f0w(nbpl)
c
      integer,parameter :: nbwl1=4
      integer :: nbpw(nbwl1)=(/2,3,5,8/)	! in nbu
c
      real(4) :: cfai(3)=(/ 1.0000,-0.3399,-0.6356/)	! FAI LCC
      real(4) :: fai1,lci,rw7
c
      real(4) :: rslope = 0.05/50000.
      real(4) :: roffset=-0.01
      real(4) :: pslope(9) =(/0.005,0.002,0.0001,0.0001,0.0001,
     &                        0.001,0.001,0.0002,1.0000/)
      real(4) :: poffset(9)=(/0.   ,0.   ,0.    ,0.    ,0.    ,
     &                        0.   ,0.   ,-0.010,0./)	! 7:TAUA_670
c      real(4) :: poffset(9)=(/0.   ,0.   ,0.    ,0.    ,0.    ,
c     &                        0.   ,-0.15,-0.010,0./)	! 7:Angstrom_670_865
c
      integer :: mnl(12)=(/0,31,59,90,120,151,181,212,243,273,304,334/)
      integer :: mnc(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: iyl( 4)=(/0,366,731,1096/)
c
      real(4),parameter :: ozn0=343.79	! DU
      real(4),parameter :: tpw0=14.186	! mm
      real(4),parameter :: epr0=1.0		! atom
      real(4):: tpw1,epr1,ozn1
      real*4 :: kwv1,kwv2,kwv3,twv0,twv1
      real*4 :: koz1,toz0,toz1
      real*4 :: kox1,kox2,kox3,tox0,tox1
c
      real(4) :: kwv(3,nbl0)=reshape( (/
     & 1.4909e-06, 0.0000e+00, 0.00, 9.8080e-07, 0.0000e+00, 0.00,
     & 3.1745e-05, 0.0000e+00, 0.00, 1.0449e-05, 0.0000e+00, 0.00,
     & 1.6566e-05, 0.0000e+00, 0.00, 1.2198e-04, 0.0000e+00, 0.00,
     & 5.6657e-05, 0.0000e+00, 0.00, 5.5299e-05, 0.0000e+00, 0.00,
     & 1.6525e-06, 0.0000e+00, 0.00, 8.0907e-05, 0.0000e+00, 0.00,
     & 7.5751e-05, 0.0000e+00, 0.00, 4.0931e-05, 0.0000e+00, 0.00,
     &-2.5585e-01, 1.4719e+00,-0.32, 4.3975e-03,-2.8813e-03, 0.06,
     & 1.7986e-02,-1.1061e-02, 0.08/),(/3,nbl0/) )
      real(4) :: koz(3,nbl0)=reshape( (/
     & 8.2534e-09, 0.0000e+00, 0.00, 2.5426e-07, 0.0000e+00, 0.00,
     & 3.0227e-06, 0.0000e+00, 0.00, 2.0641e-05, 0.0000e+00, 0.00,
     & 6.5554e-05, 0.0000e+00, 0.00, 1.1461e-04, 0.0000e+00, 0.00,
     & 4.2756e-05, 0.0000e+00, 0.00, 4.2661e-05, 0.0000e+00, 0.00,
     & 6.6933e-06, 0.0000e+00, 0.00, 1.9163e-06, 0.0000e+00, 0.00,
     & 1.8778e-06, 0.0000e+00, 0.00, 8.0493e-08, 0.0000e+00, 0.00,
     & 3.5094e-09, 0.0000e+00, 0.00, 0.0000e+00, 0.0000e+00, 0.00,
     & 0.0000e+00, 0.0000e+00, 0.00/),(/3,nbl0/) )
      real(4) :: kox(3,nbl0)=reshape( (/
     & 1.6250e-03, 0.0000e+00, 0.00, 4.2290e-05, 0.0000e+00, 0.00,
     & 4.7086e-04, 0.0000e+00, 0.00, 2.0450e-04, 0.0000e+00, 0.00,
     & 1.1597e-03, 0.0000e+00, 0.00, 5.5712e-03, 0.0000e+00, 0.00,
     & 1.9591e-03,-5.1393e-04, 1.00, 2.0069e-03,-5.2653e-04, 1.00,
     & 3.1756e-02, 2.8606e-01,-0.60, 4.4504e-05, 0.0000e+00, 0.00,
     & 4.5281e-05, 0.0000e+00, 0.00, 8.5098e-03, 0.0000e+00, 0.00,
     &-7.1343e-04, 9.0805e-04, 1.00, 4.9989e-07, 0.0000e+00, 0.00,
     & 6.5526e-08, 0.0000e+00, 0.00/),(/3,nbl0/) )
c
      integer	:: if_L2	! 0: AC, 1: add PAR, 2: LTOA/LCLR input
      real(4)	:: slope,offst
      real(4)	:: wm(0:1),wn(0:1)
      integer(4):: i,i1,i2,i3,i4
      real(4)   :: wlmin,wlmax
      integer(4):: num,iargc,nb,nb0,nb1,nbc
      integer(4):: iyr,imn,idy,ihr,imi,ise,jdy,iy1,iy2,iad
      real(4)   :: jd1
      real(4)   :: rse,rhr
      integer(4):: nlines,npixels,nl,ml,nml,m,n,im,in,iw
      integer(4):: igp,mlg,nlg,nmlg,nm,nm0
      integer(4):: md,mm,nd,nn,nmd
      integer(4):: npix0,nlin0,npix1,nlin1,idd
      real(4)   :: dgdr,dgdg,rm,rn
      integer(4):: isat,idat,nch1,nch2
      integer	:: ifsdw
      real(4)   :: avg,dev,cdv,rmx
      integer(4):: nla,mla,nal,nlo,mlo,nol
      real(4)   :: rlon,rlat,dda,ddo,rlon0,rlon1
      integer	:: ifanc,ifozn
      real(4)   :: elv1,atp1,win1
      integer(4):: ilnd
      real(4)   :: csoz,csaz,amss,ah
      real(4)   :: rd,rt,rw1,par1
      real(4)   :: sgr(nbul)
      real(4)   :: sgr0,sgr1,rsgr
      integer(4):: mo,np,npl,il,ifs
      character :: inlst(501)*48
      character(len=24) :: stt,ent,ctt
      real*4	:: area(8)		! 2:lat,lon * ul,ur,ll,lr
      integer	:: n_in,n_val,n_out
      real*4	:: r_out
      character(len=4)	:: qidx	! Good/Fair/Poor
      integer	:: ninlst,ipro0,ipro,iprol
      real*4	:: tcpu0,tcpu1	! for cpu_time
      integer*2	:: ibm			! bit mask for L3
      integer	:: norb
      character(len=1) :: res
      integer	:: vv,hh
c     ------------------------------------------------- Hokkaido & LIS
      character(len=1) :: mod44w(2801,2401)
      integer	:: ierr
c
c     ------------------------------------------------- LUT
      integer(4),parameter :: ntal=6
      integer(4),parameter :: ng1=19,ng2=19,ng3=19
      integer(4),parameter :: ngx=ng1*ng2*ng3
      integer(4),parameter :: ngd=ngx+7+ng2+1
      real(4) :: tb0(ngd,nbul)
      real(4) :: tbl(ngd,nbul,ntal,nmdl)
      real(4) :: sazg(ng1),sozg(ng2),reag(ng3)
      real(4) :: dg1,dg2,dg3
      real(4) :: wt,wta,saz1,soz1,rea1,saz2,rea2
      real(4) :: soa0,soa1,saa1,saa10,saa11,saa20,saa21
      real(4) :: rr(nbul,2)
      REAL(4) :: ra(nbul,ntal,nmdl)
      REAL(4) :: t0(nbul,ntal,nmdl)
      REAL(4) :: t1(nbul,ntal,nmdl)
      REAL(4) :: sa(nbul,ntal,nmdl)
      REAL(4) :: tb(nbul,ntal,nmdl)
      REAL(4) :: ta(nbul,ntal,nmdl)
      real(4) :: tgn(nbul)
      real(4) :: tg(nbul)
      REAL(4) :: t00(nbul)
      REAL(4) :: t10(nbul)
      REAL(4) :: sa0(nbul)
c
c     ------------------------------------------------- date_and_time
      character(LEN= 8) :: date
      character(LEN=10) :: time
      character(LEN= 5) :: zone
      integer(4)		:: idys(8)
      character(LEN=17) :: ptime
c
c     ------------------------------------------------- Allocatable
      integer(4),allocatable ::	i4buf(:)
      real(4),	 allocatable ::	r4buf(:)
      real(4),	 allocatable ::	lat(:)
      real(4),	 allocatable ::	lon(:)
      real(4),   allocatable ::	geo(:,:) ! saz,saa,soz,soa,tim,saz,saa
      integer(2),allocatable ::	i2ref(:,:)
      real(8),	 allocatable ::	tai(:)
      integer(1),allocatable ::	lwfg(:)
      integer(2),allocatable :: rrs(:,:)
      integer(2),allocatable :: par(:)
      integer(2),allocatable :: chl(:)
      integer(2),allocatable :: iop(:,:)
      integer(2),allocatable :: aot(:,:)
      integer(2),allocatable :: i2flg(:) ! nm:output flag
      integer(2),allocatable ::	fai(:)
c
      real(4),   allocatable :: r4anc(:,:,:) ! nma,ndyl,4
      integer(2),allocatable :: i2ozn(:,:)   ! nmo,ndyl
c
c     ------------------------------------------------- BRDF
      integer,parameter :: N_W=7
      integer,parameter :: N_S=6
      integer,parameter :: N_C=6
      integer,parameter :: N_N=17
      integer,parameter :: N_A=13
      integer,parameter :: N_K=5
      real(4) :: foQ0(N_A,N_N,N_C,N_S,nbwl0)
      real(4) :: foQ1(nbwl0)
      real(4) :: rgo0(N_N,N_K)
c
c     ------------------------------------------------- AC & IOP
      real(4) :: rc(nbul,2)
      real(4) :: rc0(nbul),rc1(nbul),rc2(nbul)
      real(4) :: rs0(nbul),rs1(nbul)
      real(4) :: aot1(nbul),alp,alp1,alp2
      real(4) :: rrs1(nbul)
      real(4) :: chla,aph442,adg442,bbp442,tsm,cdom
      integer :: ifiop
c
c     ---------------------------------------------------- hdf5 I/F
      CHARACTER(LEN=200)	:: filename1	! File name
      CHARACTER(LEN=200)	:: filename2	! File name
      CHARACTER(LEN=200)	:: filename3	! File name
      CHARACTER(LEN=80)		:: grp_name		! Group name
      CHARACTER(LEN=80)		:: attr_name	! Attribute name
      CHARACTER(LEN=80)		:: dsetname		! Dataset name
      CHARACTER(LEN=999)	:: attr_char	! Attribute data
      INTEGER(HID_T)		:: file_id		! File identifier
      INTEGER(HID_T)		:: file_id1		! File identifier
      INTEGER(HID_T)		:: file_id2		! File identifier
      INTEGER(HID_T)		:: file_id3		! File identifier
      INTEGER(HID_T)		:: grp_id
      INTEGER				:: class
      INTEGER				:: error		! Error flag
      INTEGER				:: error0 		! Error flag
      REAL(4)				:: r4attr(30)
      INTEGER(4)			:: i4para
      REAL(4)				:: r4para
      CHARACTER(LEN=400)	:: c1buf
c     ----------------------------------------------------
c
      call cpu_time(tcpu0)		! measure time
      fvsn=avsn(1:1)//pvsn(1:3)
c
      if(ifg_b2ss.eq.1) then
       pslope(5)=0.005
      endif
c
      if(nmdl.eq.10) then
       write(6,*) 'Aerosol model: Shettle and Fenn 1979'
       alp1=-0.013
       alp2= 1.628
      else
       write(6,*) 'Aerosol model: AERONET climatology'
       alp1=-0.148
       alp2= 2.698
      endif
c
c     ================================== Input file name
      num=iargc()
      if(num.eq.3) then
       call getarg(1,filename1)
       call getarg(2,filename2)
       call getarg(3,filename3)
      elseif(num.eq.2) then
       call getarg(1,filename1)
       call getarg(2,filename3)
       filename2=filename1
      else
       write(0,*) '> pub_sgli_Global vfile ifile ofile(or directory)'
       call exit(1)
       stop
      endif
c
      if_L2=0
      i1=len_trim(filename1)
      i2=0
      do i=1,i1-3
       if(filename1(i:i).eq.'/') i2=i
       if((filename1(i:i+3).eq.'LTOA').or.
     &    (filename1(i:i+3).eq.'LCLR')) if_L2=2
      enddo
      if(if_L2.eq.2) then
       read(filename1(i2+22:i2+25),'(2i2.2)') vv,hh
       if(filename1(i2+16:i2+16).eq.'A') then
        iad=1
       else
        iad=-1
       endif
      else
       filename1(i2+32:i2+34)='VNR'
       if(filename2.ne.'NA') filename2(i2+32:i2+34)='IRS'
       read(filename1(i2+24:i2+25),'(i2.2)') hh
       if((hh.gt.6).and.(hh.le.18)) then
        iad=-1
       else
        iad=1
       endif
      endif
      res=filename1(i2+36:i2+36)
c
      i3=len_trim(filename3)
      if(filename3(i3-1:i3).eq.'h5') then
       if_L2=1
       i4=0
       do i=1,i3
        if(filename3(i:i).eq.'/') i4=i
       enddo
       filename3(i4+ 1:i4+26)=filename1(i2+ 1:i2+26)	! date/sceneID
       filename3(i4+29:i4+30)=filename1(i2+29:i2+30)	! SG/SN
       filename3(i4+36:i4+36)=filename1(i2+36:i2+36)	! resolution
       if(filename3(i4+36:i4+36).eq.'L') then
        filename3(i4+36:i4+36)='K'
       endif
      else
       i4=i3+1
      endif
c
c     ------------------------ par range weighting
      wlmin=400.
      do nb=1,nbpl
       nb0=nbu(nbpi(nb))
       if(wvl(nb0).lt.700) then
        if(nb.lt.nbpl) then
         nb1=nbu(nbpi(nb+1))
        else
         nb1=nb0
        endif
        wlmax=(wvl(nb0)+wvl(nb1))/2.
        if(wlmax.gt.700.) wlmax=700.
        f0w(nb)=f0n(nb0)*(wlmax-wlmin)/300.
        wlmin=wlmax
       else
        f0w(nb)=f0n(nb0)*100./300.		! back-up channels
       endif
      enddo
c
      wl(1:nbul)=wvl(nbu(1:nbul))
c
c     ------------------------------------- Initialize FORTRAN interface
      CALL h5open_f(error)
      if(error.lt.0) then
       write(0,*) 'Error in h5open_f'
       call exit(1)
       stop
      endif
c
c     -------------------------------------------------------- file open
      CALL h5eset_auto_f(0,error0)
c
      CALL h5fopen_f(filename1, H5F_ACC_RDONLY_F, file_id1, error)
      if(error.lt.0) then
       write(0,*) 'No VN file = '//filename1(1:i1)
       call exit(2)
       stop
      else
       write(6,'(a)') 'Open VN file = '//filename1(1:i1)
      endif
      l1vsn=filename1(i2+38:i2+41)
c
      CALL h5fopen_f(filename2, H5F_ACC_RDONLY_F, file_id2, error)
      if(error.lt.0) then
       write(6,'(a)') 'No SW file = '//trim(filename2)
c       call exit(2)
c       stop
      else
       write(6,'(a)') 'Open SW file = '//filename2(1:i1)
      endif
c
      if(if_L2.eq.1) then
       CALL h5fopen_f(filename3, H5F_ACC_RDWR_F, file_id3, error)
       if(error.lt.0) then
        write(6,'(a)') 'No L2 file = '//filename3(1:i3)
        write(0,*) 'h5fopen_f error: '//filename3(1:i3)
        call exit(4)
        stop
       else
        write(6,'(a)') 'Add to L2 file = '//filename3(1:i3)
       endif
      else
       c1buf=filename1(i2+1:i2+41)//'.h5'
       c1buf(27:28)='L2'
       c1buf(31:35)='_NWLR'
       c1buf(38:41)=fvsn(1:4)
       filename3=filename3(1:i3)//'/'//c1buf(1:44)
       i3=len_trim(filename3)
       if(filename3(i3-8:i3-8).eq.'L') filename3(i3-8:i3-8)='K'
      endif
c
      CALL h5eset_auto_f(1,error0)
c
c     ---------------------------------- Global_attributes (file_id1)
      grp_name='Global_attributes'
      CALL h5gopen_f(file_id1,grp_name,grp_id,error)
      CALL h5eset_auto_f(0,error0)
c
      if(if_L2.eq.2) then
       attr_name='Image_start_time'
      else
       attr_name='Scene_start_time'
      endif
      call read_dattr(grp_id,attr_name,0,
     &                class,c1buf,i4para,r4para,error)
      stt=c1buf
      if(error.lt.0) then
       write(0,*) 'No Scene_start_time'
       call exit(3)
       stop
      endif
c
      if(if_L2.eq.2) then
       attr_name='Image_end_time'
      else
       attr_name='Scene_end_time'
      endif
      call read_dattr(grp_id,attr_name,0,
     &                class,c1buf,i4para,r4para,error)
      ent=c1buf
c
      if(if_L2.eq.2) then
       ctt=stt
       norb=-999
      else
       attr_name='Scene_center_time'
       CALL read_dattr(grp_id,attr_name,0,
     &                 class,c1buf,i4para,r4para,error)
       ctt=c1buf
c
       attr_name='Total_orbit_number'
       CALL read_dattr(grp_id,attr_name,0,
     &                 class,c1buf,norb,r4para,error)
      endif
c
      CALL h5eset_auto_f(1,error0)
      CALL h5gclose_f(grp_id,error)
c
      read(ctt,'(i4,2i2,2(1x,i2),1x,f6.3)') iyr,imn,idy,ihr,imi,rse
      if((mod(iyr,4).eq.0).and.(imn.gt.2)) then
       jdy=mnl(imn)+idy+1
      else
       jdy=mnl(imn)+idy
      endif
c
      dra=(0.9856003*jdy+357.528)*d2r
      dse=1.00014-0.01671*cos(dra)-0.00014*cos(2.0*dra)
      do nb=1,nbul
       f0(nb)=f0n(nbu(nb))/dse**2
      enddo
c
      kv(1:nbul)=kv0(nbu(1:nbul))
      write(6,'(a15,20f7.4)') ' Constant kv = ',kv(1:nbul)
c
      kt(1:nbul)=1.0
      kb(1:nbul)=0.0
      jd1=-1.
      if((l1vsn(1:1).eq.'0').or.(l1vsn(1:1).eq.'1').or.
     &   (l1vsn(1:1).eq.'E')) then
       jd1=int((iyr-2000)/4)*1461+iyl(mod(iyr,4)+1)+jdy
     &    +(ihr+imi/60.+rse/3600.)/24.
     &    -(int((2018-2000)/4)*1461+iyl(mod(2018,4)+1)+1)
       kt(1:nbul)=kt0(nbu(1:nbul))*jd1+1.0
       do nb=1,nbul
        nb0=nbu(nb)
        if(res.eq.'Q') then
         kb(nb)=kbq(1,nb0)+kbq(2,nb0)*jd1+kbq(3,nb0)*jd1**2
        else
         kb(nb)=kbk(1,nb0)+kbk(2,nb0)*jd1+kbk(3,nb0)*jd1**2
        endif
       enddo
       write(6,'(a4,f7.1,a4,20f7.4)') ' kt(',jd1,') = ',kt(1:nbul)
       write(6,'(a4,f7.1,a4,20f7.4)') ' kb(',jd1,') = ',kb(1:nbul)
c
      elseif(ifluce.eq.1) then
       jd1=int((iyr-2000)/4)*1461+iyl(mod(iyr,4)+1)+jdy
     &    +(ihr+imi/60.+rse/3600.)/24.
     &    -(int((2018-2000)/4)*1461+iyl(mod(2018,4)+1)+1)
       open(11,file='../lut/LUCA_error.txt',iostat=ierr,
     & access='sequential',form='formatted',status='old')
       if(ierr.eq.0) then
        i=int(jd1)+1
        if(i.lt.1) i=1
        n=0
        do while ((ierr.eq.0).and.(n.lt.i))
         read(11,'(5x,11(1x,f6.4))',iostat=ierr) (lcre(nb),nb=1,11)
         if(ierr.eq.0) n=n+1
        enddo
        write(6,'(a22,i5)') ' LUCA_error.txt line =',n
        do nb=1,nbul
         nb0=nbu(nb)
         if(nb0.le.11) then
          kt(nb)=lcre(nb0)
         endif
        enddo
       else
        do nb=1,nbul
         nb0=nbu(nb)
         if(nb0.le.11) then
          kt(nb)=kte(nb0,1)+kte(nb0,2)*jd1+kte(nb0,3)*jd1**2
     &          +kte(nb0,4)*jd1**3
         endif
        enddo
        write(6,*) 'LUCA regidual error by kte'
       endif
       write(6,'(a4,f7.1,a4,20f7.4)') ' kt(',jd1,') = ',kt(1:nbul)
      endif
c
c     ---------------------------------- Image_data attribute (file_id1)
      grp_name='Image_data'
      CALL h5gopen_f(file_id1,grp_name,grp_id,error)
      attr_name='Number_of_lines'
      CALL read_dattr(grp_id,attr_name,0,
     &                class,c1buf,nlines,r4para,error)
      attr_name='Number_of_pixels'
      CALL read_dattr(grp_id,attr_name,0,
     &                class,c1buf,npixels,r4para,error)
      attr_name='Grid_interval'
      CALL read_dattr(grp_id,attr_name,0,
     &                class,c1buf,i4para,dgdr,error)
      CALL h5gclose_f(grp_id,error)
      ml=nlines
      nl=npixels
      nml=nl*ml
c
c     ------------------------ RVS correction
      if(ifg_rvsc.gt.0) then
       if(if_L2.eq.2) then
        npixels=5000
       else
        npixels=nl
       endif
       rpix(1:npixels,1:nbl0)=1.0
       if(jdy.lt.365./2.) then
        iy1=iyr-1
        rm=jdy/365.+0.5
       else
        iy1=iyr
        rm=jdy/365.-0.5
       endif
       iy2=iy1+1
c
       if(iy1.lt.2018) iy1=2018
       if(iy2.lt.2018) iy2=2018
c
       ierr=-1
       do while ((ierr.ne.0).and.(iy1.ge.2018))
        write(c1buf,'(a34,i4.4,a4)') 
     &              '../lut/vical_pixel_N20VRS_1007_v2_',iy1,'.txt'
        open(11,file=trim(c1buf),iostat=ierr,
     &  access='sequential',form='formatted',status='old')
        iy1=iy1-1
       enddo
       if(ierr.ne.0) then
        write(6,'(a)') 'Error No RVS file (1): '//trim(c1buf)
       else
        write(6,'(a)') 'RVS file (1): '//trim(c1buf)
        do nb=1,6
         read(11,*) (rpx1(n),n=1,201)
        enddo
        close(11)
c
        ierr=-1
        do while ((ierr.ne.0).and.(iy2.ge.2018))
         write(c1buf,'(a34,i4.4,a4)') 
     &              '../lut/vical_pixel_N20VRS_1007_v2_',iy2,'.txt'
         open(12,file=trim(c1buf),iostat=ierr,
     &   access='sequential',form='formatted',status='old')
         iy2=iy2-1
        enddo
        if(ierr.ne.0) then
         write(6,'(a)') 'Error No RVS file (2): '//trim(c1buf)
        else
         write(6,'(a)') 'RVS file (2): '//trim(c1buf)
         do nb=1,6
          read(12,*) (rpx2(n),n=1,201)
          do n=1,npixels
           rp=(n-1.0)/(npixels-1.0)*200.+1.
           np=int(rp)
           rp=rp-np
           if(np.gt.200) np=200
           rpix(n,nb)=(rpx1(np)*(1.-rp)+rpx1(np+1)*rp)*(1.0-rm)
     &               +(rpx2(np)*(1.-rp)+rpx2(np+1)*rp)*rm
          enddo
c          write(6,'(21f7.4)') rm,(rpix(n,nb),n=1,npixels,250)
         enddo
         close(12)
        endif
       endif
      endif
c
c     ---------------------------------- Geometry_data attribute (file_id1)
      grp_name='Geometry_data'
      CALL h5gopen_f(file_id1,grp_name,grp_id,error)
      attr_name='Number_of_lines'
      CALL read_dattr(grp_id,attr_name,0,
     &                class,c1buf,nlines,r4para,error)
      attr_name='Number_of_pixels'
      CALL read_dattr(grp_id,attr_name,0,
     &                class,c1buf,npixels,r4para,error)
      attr_name='Grid_interval'
      CALL read_dattr(grp_id,attr_name,0,
     &                class,c1buf,i4para,dgdg,error)
c
      do i=1,8
       if    (i==1) then
        attr_name = 'Upper_left_latitude'
       elseif(i==2) then
        attr_name = 'Upper_left_longitude'
       elseif(i==3) then
        attr_name = 'Upper_right_latitude'
       elseif(i==4) then
        attr_name = 'Upper_right_longitude'
       elseif(i==5) then
        attr_name = 'Lower_left_latitude'
       elseif(i==6) then
        attr_name = 'Lower_left_longitude'
       elseif(i==7) then
        attr_name = 'Lower_right_latitude'
       elseif(i==8) then
        attr_name = 'Lower_right_longitude'
       endif
       call read_dattr(grp_id,attr_name,0, class,c1buf,i4para,area(i),
     &                  error)
      enddo
c
      CALL h5gclose_f(grp_id,error)
      mlg=nlines
      nlg=npixels
      nmlg=(mlg+1)*nlg		! +1: margin
      igp=nint(dgdg/dgdr)
c
c     --------------------------------------- Memory allocate
      allocate ( i4buf(nml) )
      allocate ( r4buf(nml) )
c
      allocate ( lat(nmlg) )
      allocate ( lon(nmlg) )
      allocate ( geo(nmlg,ngl) )
      allocate ( i2ref(nml,nbul) )
      allocate ( tai(ml) )
      allocate ( lwfg(nml) )
      allocate ( i2flg(nml) )
      i2flg(1:nml)=0
c
      allocate ( par(nml) )
      par(1:nml)=-1
c
      allocate ( rrs(nml,nbwl0) )
      rrs(1:nml,1:nbwl0)=-1
      slpw(1:nbwl0)=rslope
      offw(1:nbwl0)=roffset
c
      if(if_L2.eq.1) then
       npl=1
      else
       allocate ( chl(nml) )
       allocate ( iop(nml,3) )
       allocate ( aot(nml,2) )
       allocate ( fai(nml) )
       chl(1:nml)=-1
       iop(1:nml,1:3)=-1
       aot(1:nml,1:2)=-1
       fai(1:nml)=255
       if(if_L2.eq.2) then
        npl=nbwl0+9	! without tai
       else
        npl=nbwl0+10
       endif
      endif
c
c     -------------------------------------- Land_water_flag (file_id1)
      dsetname='/Image_data/Land_water_flag'
      CALL read2d_full(file_id1,dsetname,ml,nl, 
     &                 npix1,nlin1,i4buf,r4buf,slope,offst,error)
c
      if(error.gt.-1) then
       write(6,*) 'Read Land_water_flag'
       lwfg(1:nml)=i4buf(1:nml)
      else
       write(0,*) 'No Land_water_flag'
       lwfg(1:nml)=0
c       call exit(3)
c       stop
      endif
c
c
      if(if_L2.eq.2) then
c      ----------------------------------------- lat lon of tile data
       do mm=1,mlg
        rlat=90.-10.*vv-(mm-0.5)*dgdg
        do nn=1,nlg
         rlon=(-180.+10.*hh+(nn-0.5)*dgdg)/cos(rlat*d2r)
         nmd=(mm-1)*nlg+nn
         lat(nmd)=rlat
         lon(nmd)=rlon
        enddo
       enddo
c
      else
c      -------------------------------------- Line_tai93 (file_id1)
       dsetname='/Image_data/Line_tai93'
       CALL h5read_r8(file_id1,dsetname,ml, tai,error)
       if(error.gt.-1) then
        write(6,*) 'Read Line_tai93'
       else
        write(0,*) 'No Line_tai93'
c        call exit(3)
c        stop
       endif
c
c      ----------------------------------------- lat lon (file_id1)
       dsetname='/Geometry_data/Latitude'
       CALL read2d_full(file_id1,dsetname,mlg,nlg, 
     &                  npix0,nlin0,i4buf,lat,slope,offst,error)
c
       if(error.lt.0) then		! eqa data input
        write(0,*) 'No Latitude'
        call exit(3)
        stop
       endif
       mlg=nlin0
       nlg=npix0
c
       dsetname='/Geometry_data/Longitude'
       CALL read2d_full(file_id1,dsetname,mlg,nlg, 
     &                  npix0,nlin0,i4buf,lon,slope,offst,error)
c
       if(error.lt.0) then		! eqa data input
        write(0,*) 'No Longitude'
        call exit(3)
        stop
       endif
c     -----------------------------------------
      endif
c
c     ----------------------------------------- lwfg for Hokkaido
      if((area(1).gt. 40.).and.(area(7).lt. 46.).and.
     &   (area(4).gt.140.).and.(area(6).lt.147.)) then
       open(11,file='../lut/MOD44W_N46_E140.02801_02401_8b',
     & access='direct',recl=2801,form='unformatted',status='old',
     & iostat=ierr)
       if(ierr.eq.0) then
        write(6,*) 'Hokkaido lwfg recovered by MOD44W'
        do m=1,2401
         read(11,rec=1+m) (mod44w(n,m),n=1,2801)
        enddo
        close(11)
        do m=1,ml
         rm=(m-1.)/igp+1
         im=int(rm)
         if(im.ge.mlg) im=mlg-1
         wm(1)=rm-im
         wm(0)=1.-wm(1)
         do n=1,nl
          nm=(m-1)*nl+n
          rn=(n-1.)/igp+1
          in=int(rn)
          if(in.ge.nlg) in=nlg-1
          wn(1)=rn-in
          wn(0)=1.-wn(1)
          rlon=0.
          rlat=0.
          wta=0.
          do md=0,1
           mm=im+md
           do nd=0,1
            nn=in+nd
            nmd=(mm-1)*nlg+nn
            if((lat(nmd).ge.-90.).and.(lat(nmd).le.90.)) then
             wt=wm(md)*wn(nd)
             rlon=rlon+lon(nmd)*wt
             rlat=rlat+lat(nmd)*wt
             wta=wta+wt
            endif
           enddo
          enddo
          if(wta.gt.0) then
           rlon=rlon/wta
           rlat=rlat/wta
           if((rlat.gt. 40.).and.(rlat.lt. 46.).and.
     &        (rlon.gt.140.).and.(rlon.lt.147.)) then
            nn=nint((rlon-140.)/0.0025+1)
            mm=nint(( 46.-rlat)/0.0025+1)
            lwfg(nm)=ichar(mod44w(nn,mm))/2
           endif
          endif
         enddo
        enddo
       endif
      endif
c     ----------------------------------------- lwfg for LIS
      if((area(1).gt. 42.0).and.(area(7).lt. 41.2).and.
     &   (area(4).gt.-70.8).and.(area(6).lt.-71.7)) then
       open(11,file='../lut/MOD44W_N42_41.2_W71.7_70.8.00361_00321_8b',
     & access='direct',recl=361,form='unformatted',status='old',
     & iostat=ierr)
       if(ierr.eq.0) then
        write(6,*) 'LIS lwfg recovered by MOD44W'
        do m=1,321
         read(11,rec=1+m) (mod44w(n,m),n=1,361)
        enddo
        close(11)
        do m=1,ml
         rm=(m-1.)/igp+1
         im=int(rm)
         if(im.ge.mlg) im=mlg-1
         wm(1)=rm-im
         wm(0)=1.-wm(1)
         do n=1,nl
          nm=(m-1)*nl+n
          rn=(n-1.)/igp+1
          in=int(rn)
          if(in.ge.nlg) in=nlg-1
          wn(1)=rn-in
          wn(0)=1.-wn(1)
          rlon=0.
          rlat=0.
          wta=0.
          do md=0,1
           mm=im+md
           do nd=0,1
            nn=in+nd
            nmd=(mm-1)*nlg+nn
            if((lat(nmd).ge.-90.).and.(lat(nmd).le.90.)) then
             wt=wm(md)*wn(nd)
             rlon=rlon+lon(nmd)*wt
             rlat=rlat+lat(nmd)*wt
             wta=wta+wt
            endif
           enddo
          enddo
          if(wta.gt.0) then
           rlon=rlon/wta
           rlat=rlat/wta
           if((rlat.gt. 41.2).and.(rlat.lt. 42.0).and.
     &        (rlon.gt.-71.7).and.(rlon.lt.-70.8)) then
            nn=nint((rlon+71.7)/0.0025+1)
            mm=nint((42.0-rlat)/0.0025+1)
            lwfg(nm)=ichar(mod44w(nn,mm))/2
           endif
          endif
         enddo
        enddo
       endif
      endif
c     -----------------------------------------
c
c     ----------------------------------------- read geometry
      do nb=1,ngl
       if    (nb.eq.1) then
        file_id=file_id1
        dsetname='/Geometry_data/Sensor_zenith'
       elseif(nb.eq.2) then
        file_id=file_id1
        dsetname='/Geometry_data/Sensor_azimuth'
       elseif(nb.eq.3) then
        file_id=file_id1
        dsetname='/Geometry_data/Solar_zenith'
       elseif(nb.eq.4) then
        file_id=file_id1
        dsetname='/Geometry_data/Solar_azimuth'
       elseif(nb.eq.5) then
        file_id=file_id1
        dsetname='/Geometry_data/Obs_time'
       elseif(nb.eq.6) then
        file_id=file_id2
        dsetname='/Geometry_data/Sensor_zenith'
       elseif(nb.eq.7) then
        file_id=file_id2
        dsetname='/Geometry_data/Sensor_azimuth'
       endif
c
       CALL read2d_full(file_id,dsetname,mlg,nlg, 
     &                  npix0,nlin0,i4buf,r4buf,slope,offst,error)
       if(error.lt.0) then
        if(file_id.eq.file_id1) then
         write(0,*) '[error] No VN geometry : '//trim(dsetname)
         call exit(3)
        else
         write(6,*) '[warning] No SW geometry : '//trim(dsetname)
         i4buf(1:mlg*nlg)=-32768
        endif
       endif
c
       slpg(nb)=slope
       offg(nb)=offst
c
       do m=1,mlg
        do n=1,nlg
         nm=(m-1)*nlg+n
         if(npix0.eq.nlg) then
          nm0=nm
         elseif(npix0.lt.nlg) then
          idd=nint(1.*nlg/npix0)
          nm0=(m-1)/idd*npix0+(n-1)/idd+1
         else
          idd=nint(1.*npix0/nlg)
          nm0=(m-1)*idd*npix0+(n-1)*idd+1
         endif
         if(i4buf(nm0).lt.-32767) then
          geo(nm,nb)=-999.99
         else
          geo(nm,nb)=i4buf(nm0)*slpg(nb)+offg(nb)
         endif
        enddo
       enddo
c
      enddo		! ngl
c
c     ----------------------------------------- read reflectance
      do nb=1,nbul
       nb0=nbu(nb)
       if    (nb0.le.11) then
        file_id=file_id1
        write(dsetname,'(a17,i2.2)') '/Image_data/Lt_VN',nb0
       else
        file_id=file_id2
        write(dsetname,'(a17,i2.2)') '/Image_data/Lt_SW',nb0-11
       endif
c
       CALL read2d_full(file_id,dsetname,ml,nl, 
     &                  npix1,nlin1,i4buf,r4buf,slope,offst,error)
c
       if(error.lt.0) then
        if(file_id.eq.file_id1) then
         write(0,*) 'No reflectance : '//dsetname(1:len_trim(dsetname))
         call exit(5)
        else
         write(6,*) 'No reflectance : '//dsetname(1:len_trim(dsetname))
         i2ref(1:nml,nb)=-1
        endif
       else
        write(6,*) 'Read TOA refl.: '//dsetname(1:len_trim(dsetname))
c
        slpr(nb)=slope
        offr(nb)=offst
c
c        if(slpr(nb).lt.0.00006) then
c         isat=65535	! 1111111111111111
c        else
         isat=16383	! 0011111111111111
c         write(6,*) 'Stray-light bits are omitted'
c        endif
        do m=1,ml
         do n=1,nl
          nm=(m-1)*nl+n
          if(npix1.eq.nl) then
           nm0=nm
          elseif(npix1.lt.nl) then
           idd=nint(1.*nl/npix1)
           nm0=(m-1)/idd*npix1+(n-1)/idd+1
          else
           idd=nint(1.*npix1/nl)
           nm0=(m-1)*idd*npix1+(n-1)*idd+1
          endif
          if(i4buf(nm0).eq.65535) then	! 65535:error
           i2ref(nm,nb)=-1
          else
           idat=i4buf(nm0)
           if(nb0.eq.11) then
            if(btest(idat,15)) i2flg(nm)=ibset(i2flg(nm),7) ! straylight
           endif
           idat=iand(idat,isat)
           if(idat.ge.isat-1) then	! 16382:saturation
            i2ref(nm,nb)=-2
           else
            i2ref(nm,nb)=idat
           endif
          endif
         enddo	! n
        enddo	! m
c
       endif
      enddo		! nbul
c     -----------------------------------------------
c
      CALL h5fclose_f(file_id1, error)	! Close input file
      if(file_id2.gt.-1) CALL h5fclose_f(file_id2, error)
c
c     ---------------------------------------------- read RT tables
      call read_tbl(nmdl,nbul,nbu,sazg,sozg,reag,dg1,dg2,dg3, tb0,tbl)
c
c     ---------------------------------------------- open Ancillary(tpw,atp,win,epr)
      dda=0.5
      nla=720
      mla=361
      nal=4+1
      elv1=0.
      allocate ( r4anc(nla*mla,nal,4) )
      call read_ANC(0,iyr,jdy,nla,mla,nal, r4anc,ifanc,
     &              rlon,rlat,rhr,elv1, tpw1,atp1,win1,epr1)
c
c     ---------------------------------------------- open column ozone
      ddo=1.0
      nlo=360
      mlo=180
      nol=1+2
      allocate ( i2ozn(nlo*mlo,nol) )
      call read_TOZ(0,iyr,jdy,nol, i2ozn,ifozn, rlon,rlat,rhr, ozn1)
c     ----------------------------------------------
c
c     ----------------------------------------- rrs read
      if(if_L2.eq.1) then
       do nb=1,nbwl1
        nb1=nbpw(nb)
        if(nb1.gt.nbwl0) nb1=nbwl0
        nb0=nbu(nb1)
        dsetname='/Image_data/NWLR_'//nbw(nb1)
        CALL read2d_full(file_id3,dsetname,ml,nl, 
     &                   npix1,nlin1,i4buf,r4buf,slope,offst,error)
        if(error.lt.0) then
         write(6,*) 'No Rrs : '//trim(dsetname)
        else
         write(6,*) 'Read Rrs : '//trim(dsetname)
         slpw(nb1)=slope/f0n(nb0) ! NWLR to Rrs
         offw(nb1)=offst/f0n(nb0)
         rrs(1:nml,nb1)=i4buf(1:nml)
        endif
       enddo		! nbwl1
      endif
c     -----------------------------------------
c
      write(6,*) 'Start processing'
      n_in=0
      n_val=0
      n_out=0
      ipro0=-1
      if(nl.gt.1250) then
       iprol=20
      else
       iprol=5
      endif
c
      if(nl.gt.4000) then
       iw=iw0
      else
       iw=iw0/2
      endif
c
      do m=1,ml
c
       ipro=(m-1)*iprol/(ml-1)
       if(ipro.gt.ipro0) then
        call cpu_time(tcpu1)				! measure time
        write(6,'(a8,i2,a2,i3,f7.2,a4)') 
     &  ' Proc = ',ipro,' /',iprol,(tcpu1-tcpu0)/60.,' min'
        ipro0=ipro
       endif
c
       rm=(m-1.)/igp+1
       im=int(rm)
       if(im.ge.mlg) im=mlg-1
       wm(1)=rm-im
       wm(0)=1.-wm(1)
c
       do n=1,nl
        n_in=n_in+1
        nm=(m-1)*nl+n
        ilnd=lwfg(nm)
c
        if((ilnd.gt. 50).and.(ilnd.le.100)) then
         i2flg(nm)=ibset(i2flg(nm),1)	! land
        endif
        if((ilnd.gt.  0).and.(ilnd.lt.100)) then
         i2flg(nm)=ibset(i2flg(nm),6)	! coast
        endif
c
        nch1=0
        nch2=0
        do nb=1,nbul
         nb0=nbu(nb)
         if(i2ref(nm,nb).ne.-1) then
          nch1=nch1+1
         else
          if(nb0.le.11) nch2=nch2+1
         endif
        enddo
        if(nch1.eq.0) then
         i2flg(nm)=ibset(i2flg(nm),0)	! no observation
        endif
        if(nch2.gt.0) then
         i2flg(nm)=ibset(i2flg(nm),2)	! incomplete VNR
        endif
c
        chla  =-9.999
        aph442=-9.999
        adg442=-9.999
        bbp442=-9.999
c
c      if((m.ge.800).and.(m.le.1200).and.(n.ge.1).and.(n.le.nl)
c     &.and.(ilnd.lt.20)) then
        if((ilnd.lt.20).and.(nch1.gt.0)) then
         rn=(n-1.)/igp+1
         in=int(rn)
         if(in.ge.nlg) in=nlg-1
         wn(1)=rn-in
         wn(0)=1.-wn(1)
         wta=0.
         rlon=0.
         rlat=0.
         rhr =0.
         saz1=0.
         soz1=0.
         rea1=0.
         saa1=0.
         saz2=0.
         rea2=0.
         do md=0,1
          mm=im+md
          do nd=0,1
           nn=in+nd
           nmd=(mm-1)*nlg+nn
           if((lat(nmd).ge.-90.).and.(lat(nmd).le.90.)) then
            wt=wm(md)*wn(nd)
            rlon1=lon(nmd)
            soa1 =geo(nmd,4)
            saa11=geo(nmd,2)
            saa21=geo(nmd,7)
            if(wta.eq.0) then
             rlon0=rlon1
             soa0 =soa1
             saa10=saa11
             saa20=saa21
            else
             if(rlon1-rlon0.le.-180.) then
              rlon1=rlon1+360.
             elseif(rlon1-rlon0.gt.180.) then
              rlon1=rlon1-360.
             endif
             if(soa1-soa0.le.-180.) then
              soa1=soa1+360.
             elseif(soa1-soa0.gt.180.) then
              soa1=soa1-360.
             endif
             if(saa11-saa10.le.-180.) then
              saa11=saa11+360.
             elseif(saa11-saa10.gt.180.) then
              saa11=saa11-360.
             endif
             if(saa21-saa20.le.-180.) then
              saa21=saa21+360.
             elseif(saa21-saa20.gt.180.) then
              saa21=saa21-360.
             endif
            endif
            rlon=rlon+rlon1*wt
            rlat=rlat+lat(nmd)*wt
            rhr =rhr +geo(nmd,5)*wt
            saz1=saz1+geo(nmd,1)*wt
            soz1=soz1+geo(nmd,3)*wt
            rea1=rea1+(soa1-saa11)*wt
            saa1=saa1+saa11*wt
            saz2=saz2+geo(nmd,6)*wt
            rea2=rea2+(soa1-saa21)*wt
            wta=wta+wt
           endif
          enddo
         enddo
         if(wta.gt.0) then
          rlon=rlon/wta
          rlat=rlat/wta
          rhr =rhr /wta
          saz1=saz1/wta
          soz1=soz1/wta
          rea1=rea1/wta
          saa1=saa1/wta
          saz2=saz2/wta
          rea2=rea2/wta
          if(rea1.gt.180.) then
           rea1=rea1-360.
          elseif(rea1.le.-180.) then
           rea1=rea1+360.
          endif
          if(saa1.gt.180.) then
           saa1=saa1-360.
          elseif(saa1.le.-180.) then
           saa1=saa1+360.
          endif
          if(rea2.gt.180.) then
           rea2=rea2-360.
          elseif(rea2.le.-180.) then
           rea2=rea2+360.
          endif
          if(saz1.lt.60.) then
           n_val=n_val+1
           if(soz1.lt.89) then
c
            if(soz1.gt.75.) i2flg(nm)=ibset(i2flg(nm),11)
c
            elv1=0.
            csaz=cos(saz1*d2r)
c
            if(soz1.lt.70.) then
             csoz=cos(soz1*d2r)
            else
             ah=1.
             csoz=sin(soz1*d2r)/(6378+ah)*ah
     &           /sin(soz1*d2r-asin(6378./(6378.+ah)*sin(soz1*d2r)))
             soz1=acos(csoz)/d2r
            endif
            amss=1./csaz+1./csoz
c
            call read_ANC(1,iyr,jdy,nla,mla,nal, r4anc,ifanc,
     &                    rlon,rlat,rhr,elv1, tpw1,atp1,win1,epr1)
c
            call read_TOZ(1,iyr,jdy,nol, i2ozn,ifozn, rlon,rlat,rhr,
     &                    ozn1)
c
            call interp_tbl(nmdl,nbul,sazg,sozg,reag,epr1,tb0,tbl,
     &                      saz1,soz1,rea1, rr,ra,t0,t1,sa,ta,tb)
c
            if(win1.gt.20.) i2flg(nm)=ibset(i2flg(nm),10)	! high wind
c
            do nb=1,nbul
             nb0=nbu(nb)
             if(nb0.ge.13) then
              kwv1=kwv(1,nb0)
              kwv2=kwv(2,nb0)
              kwv3=kwv(3,nb0)
              twv0=(kwv1+kwv2*(tpw0*amss)**kwv3)*tpw0
              twv1=(kwv1+kwv2*(tpw1*amss)**kwv3)*tpw1
             else
              kwv1=kwv(1,nb0)
              twv0=kwv1*tpw0
              twv1=kwv1*tpw1
             endif
c
             if(nb0.eq.9) then
              kox1=kox(1,nb0)
              kox2=kox(2,nb0)
              kox3=kox(3,nb0)
              tox0=(kox1+kox2*(epr0*epr0*amss)**kox3)*epr0*epr0
              tox1=(kox1+kox2*(epr1*epr1*amss)**kox3)*epr1*epr1
             else
              kox1=kox(1,nb0)
              tox0=kox1*epr0*epr0
              tox1=kox1*epr1*epr1
             endif
c
             koz1=koz(1,nb0)
             toz0=koz1*ozn0
             toz1=koz1*ozn1
c
             tgn(nb)=exp(-(toz1-toz0+tox1-tox0+twv1-twv0))
             tg(nb)=tgn(nb)**amss
            enddo
c      write(6,'(20f9.1)') rhr,tpw0,tpw1,epr0,epr1,ozn0,ozn1
c      write(6,'(20f8.3)') tg(1:nbul)
c
            ifsdw=0
            do nb=1,nbul
             if(i2ref(nm,nb).ge.0) then
              if(jd1.ge.0.) then
               rd=(i2ref(nm,nb)*slpr(nb)+offr(nb))*f0(nb)/pi
               rt=(rd+kb(nb))/kt(nb)/kv(nb)/f0(nb)*pi/csoz/tg(nb)
              else
               rt=(i2ref(nm,nb)*slpr(nb)+offr(nb))/kv(nb)/csoz/tg(nb)
              endif
c
              nb0=nbu(nb)
              if((ifg_rvsc.gt.0).and.(nb0.le.6)) then	! RVS correction
               if(if_L2.eq.2) then	! LTOA
                if(saa1*iad.gt.0.) then
                 alp=asin(6378./(6378.+800.)*sin(pi-saz1*d2r))
                else
                 alp=asin(6378./(6378.+800.)*sin(pi+saz1*d2r))
                endif
                nn=nint(5000*(alp/0.754+1.)/2.)+1
                if(nn.lt.1) then
                 nn=1
                elseif(nn.gt.5000) then
                 nn=5000
                endif
                rt=rt/rpix(nn,nb0)
               else
                rt=rt/rpix(n,nb0)
               endif
              endif
c
              rc(nb,1)=rt-rr(nb,1)
              if(rc(nb,1).lt.rcmin) rc(nb,1)=rcmin
              rc(nb,2)=rt-rr(nb,2)
              if(rc(nb,2).lt.rcmin) then
               rc(nb,2)=rcmin
               if(wl(nb).lt.700.) ifsdw=ifsdw+1
              endif
             else
              rc(nb,1)=-9.9999
              rc(nb,2)=-9.9999
             endif
            enddo
c
            if(ifsdw.gt.0) i2flg(nm)=ibset(i2flg(nm),5) ! Dark pixel
c
c           ------------------------------ Sunglint
            call SUNGLINT(soz1,saz1,rea1,win1,nbul,nbu, sgr)
c
            mo=nmdl
            rsgr=1.
            do nb=1,nbul
             sgr0=sgr(nb)
             if((rc(nb,2).ge.rcmin).and.(sgr0.ge.rcmin)) then
              if(rc(nb,2).gt.0.003) then
               sgr1=(rc(nb,2)-0.003)/exp(-tb(nb,2,mo)*amss)/tg(nb)
               if(sgr1.lt.sgr0) rsgr=sgr1/sgr0
              else
               rsgr=0.
              endif
             endif
            enddo
            sgr(1:nbul)=sgr(1:nbul)*rsgr
c
            sgr0=sgr(6)
            if(sgr0.gt.tsg1) i2flg(nm)=ibset(i2flg(nm),8)	! >tsg1
            if(sgr0.gt.tsg2) i2flg(nm)=ibset(i2flg(nm),9)	! >tsg2
c
            if(sgr0.lt.tsg1) then
             n_out=n_out+1
c
c            ------------------------------ Sunglint correction for PAR
             do nb=1,nbul
              if(rc(nb,1).ge.rcmin) then
               rc0(nb)=rc(nb,1)-sgr(nb)*exp(-tb(nb,2,6)*amss)*tg(nb)
              else
               rc0(nb)=-9.999
              endif
             enddo
c
c            ------------------------------------------ Atmospheric correction
             if((if_L2.ne.1).and.(nch2.eq.0).and.(soz1.le.75.)) then
c
c             ------------------------------ Sunglint correction
              do nb=1,nbul
               if(rc(nb,2).ge.rcmin) then
                if(ifg_sung.eq.1) then
                 rc1(nb)=rc(nb,2)-sgr(nb)*exp(-tb(nb,3,mo)*amss)*tg(nb)
                else
                 rc1(nb)=rc(nb,2)
                endif
                rc2(nb)=rc1(nb)
     &                 /(t0(nb,1,1)*t1(nb,1,1)+sa(nb,1,1)*rc1(nb))
               else
                rc1(nb)=-9.999
                rc2(nb)=-9.999
               endif
              enddo
c
c             ---------------------------- rt deviation of nb=3
              num=0
              avg=0.
              dev=0.
              rmx=0.
              nb=3
              do md=-iw,iw
               mm=m+md
               if((mm.gt.0).and.(mm.le.ml)) then
                do nd=-iw,iw
                 nn=n+nd
                 if((nn.gt.0).and.(nn.le.nl)) then
                  nmd=(mm-1)*nl+nn
                  if((i2ref(nmd,nb).ge.0).and.(lwfg(nmd).lt.20)) then
                   rt=i2ref(nmd,nb)*slpr(nb)+offr(nb)
                   num=num+1
                   avg=avg+rt
                   dev=dev+rt*rt
                   if(rt.gt.rmx) rmx=rt
                  endif
                 endif
                enddo
               endif
              enddo
              if(num.gt.iw*2) then
               avg=avg/num
               rt=i2ref(nm,nb)*slpr(nb)+offr(nb)
               cdv=(dev+rt*rt)/(num+1)-avg**2
               if(cdv.gt.0) then
                cdv=sqrt(cdv)/csoz/tg(nb)
               else
                cdv=0.
               endif
               dev=dev/num-avg**2
               if(dev.gt.0) then
                dev=sqrt(dev)/csoz/tg(nb)
               else
                dev=0.
               endif
              else
               dev=0.
               cdv=0.
              endif
              rmx=rmx/kv(nb)/csoz/tg(nb)-rr(nb,2)
              if((rmx.ge.rcmin).and.(ifg_sung.eq.1)) then
               rmx=rmx-sgr(nb)*exp(-tb(nb,4,6)*amss)*tg(nb)
              endif
c
c             ------------------------------------------ Cloud screening
c      write(6,'(9f8.4)') rc1(3),rc1(4),cdv,rc1(7),rc1(10),rmx,dev
              if((rc1(3).lt.0.20).and.(rc1(4).lt.0.20).and.
     &           (ifsdw.eq.0).and.(cdv.lt.0.03).and.
     &           (rc1(7).gt.0.).and.(rc1(10).gt.0.).and.
     &           ((rc1(7).lt.0.09).or.(rc1(10).lt.0.09))) then
c
               if((rmx.ge.0.20).or.
     &            ((rc1(10).ge.0.03).and.(dev.ge.0.005)).or.
     &            (abs(rc1(7)-rc1(8)).gt.sgr(7)+0.03)) then
                i2flg(nm)=ibset(i2flg(nm),4)	! near cloud
               endif
c
c              -------------------- First guess
               call rw_sim_lci(nbul,rc2, rs1)
c
c              -------------------- FAI
               if((rc2(10).ge.0.).and.(rc2(13).ge.0.).and.
     &            (rc2( 7).ge.0.).and.
     &            (rc2( 4).ge.0.).and.(rc2( 6).ge.0.)) then
                lci=rc2(4)-1.2226*rc2(6)+0.1929*rc2(13)		! VN4,VN6,SW03
                rw7=0.00252-0.66871*lci-12.90886*lci**2-375.50607*lci**3
                if(rw7.lt.0.003) then
                 rw7=0.003
                elseif(rw7.gt.rc2(7)) then
                 rw7=rc2(7)
                endif
                rw7=rw7-0.003
c
                fai1=cfai(1)*rc2(10)
     &              +cfai(2)*rc2(13)
     &              +cfai(3)*(rc2(7)-rw7)
c                -------------------------- comment out*
c                if(fai1.gt.rs1(10)) then
c                 rs1(10)=fai1
c                 rs1(11)=fai1
c                endif
c                if(rw7.gt.rs1(7)) then
c                 rs1(7)=rw7
c                 rs1(8)=rw7
c                endif
c                -------------------------- 
               else
                fai1=-9.999
               endif
c
c              --------------------------------------- AOT & rw
               call estimate_aot(nmdl,nbul,nbu,wl,ra,t0,t1,sa,ta,tb,tg,
     &             rc1,rs1,sgr0,amss,nbal,nba,ifngc, aot1,rs0)
c              ---------------------------------------
c
               if((aot1(10).gt.0.05).and.(aot1(7).gt.0.05)) then
                alp=-log(aot1(7)/aot1(10))/log(wl(7)/wl(10))
                if((alp.lt.alp1).or.(alp.gt.alp2)) then
                 i2flg(nm)=ibset(i2flg(nm),13)		! out of models
                endif
               else
                alp=-9.999
                if(aot1(10).lt.-9.) i2flg(nm)=ibset(i2flg(nm),13)
               endif
c
c              -------------------- BRDF corr & IOP
               rrs1(1:nbwl0)=rs0(1:nbwl0)/pi	! Rs -> Rrs
               ifiop=0
               call cal_iop(ifiop,rrs1, chla,aph442,adg442,bbp442)
c
               if((ifg_brdf.gt.0).and.(chla.gt.0.001)) then
                call brdf_oc(ifg_brdf,nbwl0,nbul,wl,saz1,soz1,rea1,chla,
     &                       win1, rgo0,foQ0, foQ1)
                rrs1(1:nbwl0)=rrs1(1:nbwl0)/foQ1(1:nbwl0)
                ifiop=1
                call cal_iop(ifiop,rrs1, chla,aph442,adg442,bbp442)
               endif
c      write(6,'(a5,99f8.4)') 'rrs1=',rrs1(1:nbwl0)
c      write(6,'(99f8.3)') chla,aph442,adg442,bbp442
c              --------------------
c
               if((rs0(7).gt.0.008).or.(chla.gt.60.))
     &                             i2flg(nm)=ibset(i2flg(nm),15)	! Case2
c
               if(rs0(6).lt.0.002) i2flg(nm)=ibset(i2flg(nm),5)	! Dark pixel
c
c              ------------------------------------------ Set outputs
c
c              -------------------- Rrs
               nbc=0
               do nb=1,nbwl0
                if(rrs1(nb).lt.0) nbc=nbc+1
                idat=nint((rrs1(nb)-roffset)/rslope)
                if(idat.lt.0) then
                 idat=-1
                elseif(idat.gt.32767) then
                 if(idat.gt.65534) idat=65534
                 idat=idat-65536
                endif
                rrs(nm,nb)=idat
               enddo
c
               if(nbc.gt.0) i2flg(nm)=ibset(i2flg(nm),14) ! negative Rrs
c
c              -------------------- CHLA
               idat=nint((chla-poffset(2))/pslope(2))
               if(idat.lt.0) then
                idat=-1
               elseif(idat.gt.32767) then
                if(idat.gt.65534) idat=65534
                idat=idat-65536
               endif
               chl(nm)=idat
c
c              -------------------- IOP
               idat=nint((aph442-poffset(3))/pslope(3))
               if(idat.lt.0) then
                idat=-1
               elseif(idat.gt.32767) then
                if(idat.gt.65534) idat=65534
                idat=idat-65536
               endif
               iop(nm,1)=idat
c
               if(ifg_b2ss.eq.0) then
                idat=nint((adg442-poffset(4))/pslope(4))
               else
                if(adg442.gt.0) then
                 cdom=adg442*(1.54332/0.98481)*(1.-0.1)	! 442->412, adg->ag
                 idat=nint((cdom-poffset(4))/pslope(4))
                else
                 idat=-1
                endif
               endif
               if(idat.lt.0) then
                idat=-1
               elseif(idat.gt.32767) then
                if(idat.gt.65534) idat=65534
                idat=idat-65536
               endif
               iop(nm,2)=idat
c
c              -------------------- TSM
               if(ifg_b2ss.eq.0) then
                idat=nint((bbp442-poffset(5))/pslope(5))
               else
                if(bbp442.gt.0) then
c                 tsm=68.631*bbp442+1811.666*bbp442**2
                 tsm=64.985*bbp442+1150.731*bbp442**2	! Valente 2019 (bbp442>0.005)
                 idat=nint((tsm-poffset(5))/pslope(5))
                else
                 idat=-1
                endif
               endif
               if(idat.lt.0) then
                idat=-1
               elseif(idat.gt.32767) then
                if(idat.gt.65534) idat=65534
                idat=idat-65536
               endif
               iop(nm,3)=idat
c
c              -------------------- AOT
               idat=nint((aot1(10)-poffset(6))/pslope(6))
               if(idat.lt.0) then
                idat=-1
               elseif(idat.gt.32767) then
                if(idat.gt.65534) idat=65534
                idat=idat-65536
               endif
               aot(nm,1)=idat
c
               if(aot1(10).gt.0.5) i2flg(nm)=ibset(i2flg(nm),12) ! AOT
c
               if(poffset(7).lt.0.) then
                idat=nint((alp-poffset(7))/pslope(7))		! Angstr_670_865
               else
                idat=nint((aot1(7)-poffset(7))/pslope(7))	! TAUA_670
               endif
               if(idat.lt.0) then
                idat=-1
               elseif(idat.gt.32767) then
                if(idat.gt.65534) idat=65534
                idat=idat-65536
               endif
               aot(nm,2)=idat
c
c              -------------------- FAI
               if(fai1.gt.-9.) then
                if    ( btest(i2flg(nm),4) ) then	! mask near cloud
                 idat=255
                elseif( btest(i2flg(nm),6) ) then	! mask coast
                 idat=255
                else
                 idat=nint((fai1-poffset(8))/pslope(8))
                 if(idat.lt.0) then
                  idat=0
                 elseif(idat.gt.254) then
                  idat=254
                 endif
                endif
                fai(nm)=idat
               else
                fai(nm)=255
               endif
c              ------------------------------------------
c
              else
               i2flg(nm)=ibset(i2flg(nm),3)	! cloud or ice
              endif
c
             endif		! ocean color processing
c
c            ------------------------------------------ PAR
             do nb=1,nbul
              rw1=-7.4095E-03+1.2507E+09*wl(nb)**(-4)
              if(rw1.gt.0.035) rw1=0.035
              if(rw1.lt.rcmin) rw1=rcmin
              rs0(nb)=rw1
             enddo
c
             do nb=1,nbpl
              nb0=nbpi(nb)
              if(nb0.le.8) then
               nb1=nb0
               if(nb1.gt.nbwl0) nb1=nbwl0
               if(rrs(nm,nb1).ge.0) then
                rs0(nb0)=(rrs(nm,nb1)*slpw(nb1)+offw(nb1))*pi
                if(rs0(nb0).lt.rcmin) rs0(nb0)=rcmin
               elseif(rrs(nm,nb1).lt.-2) then
                rs0(nb0)=((rrs(nm,nb1)+65536)*slpw(nb1)+offw(nb1))*pi
                if(rs0(nb0).lt.rcmin) rs0(nb0)=rcmin
               endif
              endif
             enddo
c
             t00(1:nbul)=t0(1:nbul,1,1)
             t10(1:nbul)=t1(1:nbul,1,1)
             sa0(1:nbul)=sa(1:nbul,1,1)
c
             do nb=1,nbul
              if(rr(nb,2).gt.rr(nb,1)) then	! put surface reflection
               rs1(nb)=rs0(nb)+(rr(nb,2)-rr(nb,1))/t00(nb)/t10(nb)
              else
               rs1(nb)=rs0(nb)
              endif
             enddo
c
             call cal_par(nbul,nbpl,nbpi,iyr,jdy,rhr,rlat,csoz,rc0,rs1,
     &                    tgn,t00,t10,sa0,f0w, par1)
             if(par1.ge.0) then
              par(nm)=nint((par1-poffset(1))/pslope(1))
             endif
c            ------------------------------------------
c
c      write(6,'(i4,99f8.3)') n,par1
c      write(6,'(99f8.4)') rc(1:nbul,1)
c      write(6,'(99f8.4)') rc0(1:nbul)
c      write(6,'(99f8.4)') rs1(1:nbul)
c
            endif	! sgr0<tsg1
           else
            par(nm)=0
           endif	! soz
          endif		! saz1<60
         endif		! wta
        endif		! ocean
       enddo	! nl
      enddo		! ml
c
c     ---------------------------------------- near cloud (by flag-3)
c      if(if_L2.ne.1) then
c       write(6,*) 'Near cloud (by flag-3)'
c       do m=1,ml
c        do n=1,nl
c         nm=(m-1)*nl+n
c         if((btest(i2flg(nm),1).eqv..false.).and.
c     &      (btest(i2flg(nm),2).eqv..false.).and.
c     &      (btest(i2flg(nm),3).eqv..false.)) then
c          i=0
c          do md=-iw,iw
c           mm=m+md
c           if((mm.gt.0).and.(mm.le.ml)) then
c            do nd=-iw,iw
c             nn=n+nd
c             if((nn.gt.0).and.(nn.le.nl)) then
c              nmd=(mm-1)*nl+nn
c              if( btest(i2flg(nmd),3) ) i=1
c             endif
c            enddo
c           endif
c          enddo
c          if(i.gt.0) i2flg(nm)=ibset(i2flg(nm),4)
c         endif
c        enddo
c       enddo
c      endif
c     ----------------------------------------
c
c     ------------------------------------------------------- Attributes
      if(if_L2.ne.1) then
c
       CALL h5fcreate_f(filename3, H5F_ACC_TRUNC_F, file_id3, error)
       if(error.lt.0) then
        write(0,*) 'Cannot make L2 file = '//filename3(1:i3)
        call exit(4)
        stop
       endif
       write(6,*) 'Create L2 file = '//filename3(1:i3)
c
       c1buf=filename3(i4+1:i3)
       if(if_L2.eq.2) then
        il=12
       else
        il=16
       endif
       call gl_attr(file_id3,il,c1buf,stt,ent,ctt,avsn,pvsn,norb,error)
c
c      ------------------------------------ Product quality index
       if(n_in.eq.0) then
        qidx='Fair'
        write(6,'(a8,i10)') ' n_in  =',n_in
       elseif(n_val.eq.0) then
        qidx='Fair'
        write(6,'(a8,i10)') ' n_val =',n_val
       else
        r_out=1.*n_out/nl/ml
        write(6,'(a8,f8.5)') ' r_out =',r_out
        if(r_out.lt.0.0001) then
         qidx='Fair'
        else
         qidx='Good'
        endif
       endif
c      ------------------------------------ get processing date
       write(6,*) 'Get processing date'
       call date_and_time(date,time,zone,idys)
       rhr=idys(5)+idys(6)/60.+idys(7)/3600.+idys(8)/3600000.
     &    -idys(4)/60.
       if(rhr.lt.0) then
        rhr=rhr+24.
        idys(3)=idys(3)-1
        if(idys(3).lt.1) then
         idys(2)=idys(2)-1
         if(idys(2).lt.1) then
          idys(3)=31
          idys(2)=12
          idys(1)=idys(1)-1
         elseif((mod(idys(1),4).eq.0).and.(idys(2).eq.2)) then
          idys(3)=29
         else
          idys(3)=mnc(idys(2))
         endif
        endif
       endif
       ihr=int(rhr)
       imi=int((rhr-ihr)*60)
       ise=int(((rhr-ihr)*60-imi)*60)
       write(ptime,'(i4,2i2.2,3(a1,i2.2))')
     & idys(1),idys(2),idys(3),' ',ihr,':',imi,':',ise
c      write(6,'(a16,a17)') 'Processing UT = ',ptime
c      ------------------------------------ input list
       ninlst=2
       inlst(1)=filename1(i2+1:len_trim(filename1))
       inlst(2)=filename1(i2+1:len_trim(filename2))
c
       call sp_attr(file_id3,nl,ml,dgdr,nlg,mlg,dgdg,ptime,area,
     &              ninlst,inlst,qidx,ifanc,ifozn, error)
c
      else
       write(6,*) 'Write to '//filename3(1:i3)
      endif
c
c     ------------------------------------------------------- datasets
      do np=1,npl	! rrs,PAR,Chl,aph,adg,bbp,aot8,alp,fai,QA_flag,tai
       ibm=0
       ibm=ibset(ibm,0)	! data lack
       ibm=ibset(ibm,1)	! land
       ibm=ibset(ibm,2)	! incomplete band
       if(np.eq.1) then
c       ------------------------------------- PAR
        dsetname='/Image_data/PAR'
        attr_char='Photosynthetically Available Radiation'
        r4attr(3)=pslope(1)
        r4attr(4)=poffset(1)
        il=10
        ifs=2		! uint16
        do nm=1,nml
         if(par(nm).lt.0) then
          i4buf(nm)=par(nm)+65536
         else
          i4buf(nm)=par(nm)
         endif
        enddo
        r4attr(8)=dgdr
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.le.1+nbwl0) then
c       ------------------------------------- NWLR
        nb=np-1
        dsetname='/Image_data/NWLR_'//nbw(nb)
        attr_char='Normalized Water Leaving Radiance(NWLR) at '//nbw(nb)
     &//' nm [W/m^2/sr/um] = DN * Slope + Offset; Remote Sensing Reflect
     &ance(Rrs) at '//nbw(nb)//' [sr^-1] = DN * Rrs_slope + Rrs_offset'
        r4attr(3)=rslope*f0n(nbu(nb))
        r4attr(4)=roffset*f0n(nbu(nb))
        do nm=1,nml
         if(rrs(nm,nb).lt.0) then
          i4buf(nm)=rrs(nm,nb)+65536
         else
          i4buf(nm)=rrs(nm,nb)
         endif
        enddo
        r4attr(8)=dgdr
        r4attr(11)=wl(nb)
        r4attr(13)=bw(nbu(nb))
        r4attr(15)=f0(nb)
        r4attr(17)=rslope
        r4attr(18)=roffset
        ifs=2		! uint16
        il=19
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+2) then
c       ------------------------------------- Chl
        dsetname='/Image_data/CHLA'
        attr_char='Chlorophyll-a concentration'
        r4attr(3)=pslope(np-nbwl0)
        r4attr(4)=poffset(np-nbwl0)
        il=10
        ifs=2		! uint16
        do nm=1,nml
         if(chl(nm).lt.0) then
          i4buf(nm)=chl(nm)+65536
         else
          i4buf(nm)=chl(nm)
         endif
        enddo
        r4attr(8)=dgdr
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+3) then
c       ------------------------------------- aph
        dsetname='/Image_data/aph_442'
        attr_char='Absorption coefficient of phytoplankton at 442 nm'
        r4attr(3)=pslope(np-nbwl0)
        r4attr(4)=poffset(np-nbwl0)
        il=10
        ifs=2		! uint16
        do nm=1,nml
         if(iop(nm,1).lt.0) then
          i4buf(nm)=iop(nm,1)+65536
         else
          i4buf(nm)=iop(nm,1)
         endif
        enddo
        r4attr(8)=dgdr
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+4) then
c       ------------------------------------- adg
        if(ifg_b2ss.eq.0) then
         dsetname='/Image_data/adg_442'
         attr_char='Absorption coefficient of detritus+CDOM at 442 nm'
        else
         dsetname='/Image_data/CDOM'
         attr_char='Colored dissolved organic matter (CDOM) at 412nm = D
     &N * Slope + Offset [m^-1]'
        endif
        r4attr(3)=pslope(np-nbwl0)
        r4attr(4)=poffset(np-nbwl0)
        il=10
        ifs=2		! uint16
        do nm=1,nml
         if(iop(nm,2).lt.0) then
          i4buf(nm)=iop(nm,2)+65536
         else
          i4buf(nm)=iop(nm,2)
         endif
        enddo
        r4attr(8)=dgdr
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+5) then
c       ------------------------------------- bbp or TSM
        if(ifg_b2ss.eq.0) then
         dsetname='/Image_data/bbp_442'
         attr_char='Backscattering coefficient of particles at 442 nm'
        else
         dsetname='/Image_data/TSM'
         attr_char='Total suspended matter (TSM) = DN * Slope + Offset [
     &g m^-3]'
        endif
        r4attr(3)=pslope(np-nbwl0)
        r4attr(4)=poffset(np-nbwl0)
        il=10
        ifs=2		! uint16
        do nm=1,nml
         if(iop(nm,3).lt.0) then
          i4buf(nm)=iop(nm,3)+65536
         else
          i4buf(nm)=iop(nm,3)
         endif
        enddo
        r4attr(8)=dgdr
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+6) then
c       ------------------------------------- TAUA_865
        dsetname='/Image_data/TAUA_865'
        attr_char='Aerosol optical thickness at 865nm'
        r4attr(3)=pslope(np-nbwl0)
        r4attr(4)=poffset(np-nbwl0)
        il=10
        ifs=2		! uint16
        do nm=1,nml
         if(aot(nm,1).lt.0) then
          i4buf(nm)=aot(nm,1)+65536
         else
          i4buf(nm)=aot(nm,1)
         endif
        enddo
        r4attr(8)=dgdr
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+7) then
c       ------------------------------------- TAUA_670 or Angstrom_670_865
        r4attr(3)=pslope(np-nbwl0)
        r4attr(4)=poffset(np-nbwl0)
        if(r4attr(4).lt.0.) then
         dsetname='/Image_data/Angstrom_670_865'
         attr_char='Aerosol angstrom exponent between 670 and 865 nm'
         ibm=ibset(ibm,13)	! out of models
        else
         dsetname='/Image_data/TAUA_670'
         attr_char='Aerosol optical thickness at 670nm'
        endif
        il=10
        ifs=2		! uint16
        do nm=1,nml
         if(aot(nm,2).lt.0) then
          i4buf(nm)=aot(nm,2)+65536
         else
          i4buf(nm)=aot(nm,2)
         endif
        enddo
        r4attr(8)=dgdr
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+8) then
c       ------------------------------------- FAI
        dsetname='/Image_data/FAI'
        attr_char='Floating algae index'
        r4attr(3)=pslope(np-nbwl0)
        r4attr(4)=poffset(np-nbwl0)
        il=10
        ifs=1		! uint8
        i4buf(1:nml)=fai(1:nml)
        r4attr(8)=dgdr
        ibm=ibset(ibm,3)	! cloud or ice
        ibm=ibset(ibm,4)	! near cloud
        ibm=ibset(ibm,6)	! coast
        call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+9) then
c       ------------------------------------- QA_flag (file_id3)
        dsetname='/Image_data/QA_flag'
        write(attr_char,'(16a28)')
     &   'Bit00:no observation data  ;',
     &   'Bit01:land pixel           ;',
     &   'Bit02:incomplete VNR bands ;',
     &   'Bit03:cloud or ice         ;',
     &   'Bit04:near cloud (+-2pix)  ;',
     &   'Bit05:dark pixel           ;',
     &   'Bit06:coast pixel          ;',
     &   'Bit07:straylight flag      ;',
     &   'Bit08:sunglint mask>0.16   ;',
     &   'Bit09:sunglint flag>0.01   ;',
     &   'Bit10:wind speed>20m/s     ;',
     &   'Bit11:soz>75               ;',
     &   'Bit12:taua>0.5             ;',
     &   'Bit13:out of aerosol models;',
     &   'Bit14:negative nLw         ;',
     &   'Bit15:turbid Case-2 water  ;'
         r4attr(3)=pslope(np-nbwl0)
         r4attr(4)=poffset(np-nbwl0)
         il=9
         ifs=2			! uint16
         do nm=1,nml
          if(i2flg(nm).lt.0) then
           i4buf(nm)=i2flg(nm)+65536
          else
           i4buf(nm)=i2flg(nm)
          endif
         enddo
         call h5_write(nl,ml,file_id3,dsetname,i4buf,r4buf,ifs,
     &                 attr_char,il,r4attr,ibm, error)
c
       elseif(np.eq.nbwl0+10) then
c       -------------------------------------- Line_tai93 (->file_id3)
        dsetname='/Image_data/Line_tai93'
        CALL h5_write_r8(file_id3,dsetname,ml,tai, error)
c
       endif
c
       if(error.lt.0) then
        write(0,*) ' Error write : '//trim(dsetname)
        call exit(6)
        stop
       else
        write(6,*) 'Write '//trim(dsetname)
       endif
      enddo	! np
c
c     -----------------------------------------Geo
      if(if_L2.eq.0) then
       nmlg=nlg*mlg
       r4attr(3)=1.
       r4attr(4)=0.
       r4attr(8)=igp
       ifs=0			! uint16
       il=9
       do nb=1,2
        if    (nb.eq.1) then
         dsetname='/Geometry_data/Latitude'
         attr_char='Latitude in degrees'
         r4buf(1:nmlg)=lat(1:nmlg)
        elseif(nb.eq.2) then
         dsetname='/Geometry_data/Longitude'
         attr_char='Longitude in degrees'
         r4buf(1:nmlg)=lon(1:nmlg)
        endif
        call h5_write(nlg,mlg,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
       enddo
      endif
c
      if(if_L2.ne.1) then
       nmlg=nlg*mlg
       r4attr(8)=igp
       ifs=-2			! uint16
       il=9
       do nb=1,5
        if    (nb.eq.1) then
         dsetname='/Geometry_data/Sensor_zenith'
         attr_char='Sensor zenith angle in degrees'
        elseif(nb.eq.2) then
         dsetname='/Geometry_data/Sensor_azimuth'
         attr_char='Sensor azimuth angle in degrees'
        elseif(nb.eq.3) then
         dsetname='/Geometry_data/Solar_zenith'
         attr_char='Solar zenith angle in degrees'
        elseif(nb.eq.4) then
         dsetname='/Geometry_data/Solar_azimuth'
         attr_char='Solar azimuth angle in degrees'
        elseif(nb.eq.5) then
         dsetname='/Geometry_data/Obs_time'
         attr_char='Observation time in hours'
        endif
        do nm=1,nmlg
         idat=nint((geo(nm,nb)-offg(nb))/slpg(nb))
         if((idat.lt.-32768).or.(idat.gt.32767)) idat=-32768
         i4buf(nm)=idat
        enddo
        r4attr(3)=slpg(nb)
        r4attr(4)=offg(nb)
        call h5_write(nlg,mlg,file_id3,dsetname,i4buf,r4buf,ifs,
     &                attr_char,il,r4attr,ibm, error)
       enddo
      endif
c
c     -------------------- Memory deallocate
      write(6,*) 'Deallocate Arrays'
      deallocate ( i4buf )
      deallocate ( r4buf )
      deallocate ( rrs )
      deallocate ( i2ref )
      deallocate ( tai )
      deallocate ( lwfg )
      deallocate ( lat,lon,geo )
      deallocate ( par )
      deallocate ( r4anc,i2ozn )
      if(if_L2.ne.1) then
       deallocate ( chl,iop,aot,fai )
      endif
c     --------------------
c
      CALL h5fclose_f(file_id3, error)	! Close output file
c
      CALL h5close_f(error)				! Close FORTRAN interface
      if(error.lt.0) then
       write(0,*) 'Error in h5fclose_f(file_id3)'
       call exit(1)
       stop
      endif
c
      write(6,*) 'Finish'
      call exit(100)
c
      end
c
c=======================================================================
c
c=======================================================================
      SUBROUTINE read_tbl(nmdl,nbsl,nbu,sazg,sozg,reag,dg1,dg2,dg3,
     &                    tb0,tbl)
c
      implicit none
c
c    6:WL,F0,rsf,n_saz,n_soz,n_rea
c    ng1+ng2+ng3: geo grids
c    ngx: rtoa
c    7:Total, Particle, and Rayleigh optical thickness and 
c     Single scattering albedo,PRS,TMP,tpw
c    ng2:total transmittance at soz
c    1: pherical albedo
c
      integer,parameter :: nbsl0=17
      integer,parameter :: ntal=6
      integer,parameter :: ng1=19,ng2=19,ng3=19
      integer,parameter :: nsk=6+ng1+ng2+ng3
      integer,parameter :: ngx=ng1*ng2*ng3
      integer,parameter :: ngd=ngx+7+ng2+1
      integer,parameter :: ngi=nsk+ngd
c
      integer :: nmdl,nbsl
      integer :: nbu(nbsl)
      real(4) :: sazg(ng1),sozg(ng2),reag(ng3)
      real(4) :: dg1,dg2,dg3
      real(4) :: tb0(ngd,nbsl)
      real(4) :: tbl(ngd,nbsl,ntal,nmdl)
c
      integer :: i,j,it,nb,nb0,ierr,n1,i2,im
      integer :: itl,iml
      character :: ftbl*200
      real(4) :: r4buf(ngi)
c
      do j=1,2
       if(j.eq.1) then
c        ftbl='Rayleigh_pstar4I_SGLI0_lmb_flt_t00_5deg'
        ftbl='LUT_pstar4I_SGLI0_lmb_flt_t00_5deg_buoy'
        iml=1
        itl=1
       elseif(j.eq.2) then
        if(nmdl.gt.9) then
         ftbl='LUT_pstar4I_SGLI0_r00_t00_01_03_06_10_15_g01'
        else
c         ftbl='LUT_pstar4I_SGLI0_fl5_t00_01_03_06_10_15_o01' ! Kosa*0
c         ftbl='LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o01' ! Kosa*0.05
         ftbl='LUT_pstar4I_SGLI0_f5j_t00_01_03_06_10_15_o01_buoy'
        endif
        iml=nmdl
        itl=ntal
       endif
c       i2=len_trim(ftbl)
       i2=len_trim(ftbl)-5
c
c     ----------------------------------------------- table grid
       open(11,file='../lut/'//ftbl,iostat=ierr,
     & access='direct',form='unformatted',recl=4*ngi,status='old')
       if(ierr.ne.0) then
        open(11,file='lut/'//ftbl,iostat=ierr,
     &  access='direct',form='unformatted',recl=4*ngi,status='old')
        if(ierr.ne.0) then
         write(6,'(a)') 'No LUT = '//trim(ftbl)
         call exit(10)
        endif
       endif
       read(11,rec=1) (r4buf(i),i=1,ngi)
       close(11)
       do i=1,ng1
        sazg(i)=180.-r4buf(6+i)
       enddo
       do i=1,ng2
        sozg(i)=r4buf(6+ng1+i)
       enddo
       do i=1,ng3
        reag(i)=180.-r4buf(6+ng1+ng2+i)
        if(reag(i).gt.180) then
         reag(i)=reag(i)-360.
        endif
       enddo
       dg1=abs(sazg(2)-sazg(1))
       dg2=abs(sozg(2)-sozg(1))
       dg3=abs(reag(2)-reag(1))
c
c      ----------------------------------------------- read data
       do im=1,iml
        if(j.eq.2) then
         write(ftbl(i2-1:i2),'(i2.2)') im
        endif
        write(6,'(a)') 'read LUT = '//trim(ftbl)
        open(11,file='../lut/'//ftbl,iostat=ierr,
     &  access='direct',form='unformatted',recl=4*ngi,status='old')
        if(ierr.ne.0) then
         open(11,file='lut/'//ftbl,iostat=ierr,
     &   access='direct',form='unformatted',recl=4*ngi,status='old')
        endif
        do it=1,itl
         do nb=1,nbsl
          nb0=nbu(nb)
          n1=(it-1)*nbsl0+nb0
          read(11,rec=n1) (r4buf(i),i=1,ngi)
          do i=1,ngd
           if(j.eq.1) then
            tb0(i,nb)=r4buf(nsk+i)
           else
            tbl(i,nb,it,im)=r4buf(nsk+i)
           endif
          enddo
         enddo
        enddo
        close(11)
       enddo
c      -----------------------------------------------
      enddo
c     -----------------------------------------------
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE interp_tbl(nmdl,nbsl,sazg,sozg,reag,epr1,tb0,tbl,
     &                      saz,soz,rea, rr,ra,t0,t1,sa,ta,tb)
c
      implicit none
c
      integer,parameter :: ntal=6
      integer,parameter :: ng1=19,ng2=19,ng3=19
      integer,parameter :: ngx=ng1*ng2*ng3
      integer,parameter :: ix=ngx+7		! skip for total trans
      integer,parameter :: ngd=ix+ng2+1
      integer,parameter :: igl=1
c
      integer :: nmdl,nbsl
      real(4) :: sazg(ng1),sozg(ng2),reag(ng3)
      real(4) :: epr1
      real(4) :: tb0(ngd,nbsl)
      real(4) :: tbl(ngd,nbsl,ntal,nmdl)
      real(4) :: saz,soz,rea
      real(4) :: rr(nbsl,2)		! 2: land,ocean surface
      real(4) :: ra(nbsl,ntal,nmdl)
      real(4) :: t0(nbsl,ntal,nmdl)
      real(4) :: t1(nbsl,ntal,nmdl)
      real(4) :: sa(nbsl,ntal,nmdl)
      real(4) :: ta(nbsl,ntal,nmdl)
      real(4) :: tb(nbsl,ntal,nmdl)
c
      real(4) :: dg1,dg2,dg3
      real(4) :: rg1(-igl:igl),rg2(-igl:igl),rg3(-igl:igl)
      integer :: ig1(-igl:igl),ig2(-igl:igl),ig3(-igl:igl)
      integer :: i,im,iml,it,itl,nb,ig0,id0,id3,id2,id1
      integer :: ix1,ix2
      real(4) :: wbas,wmin,rg0,wgt,wt
      real(4) :: rt0,rt1,ra0,tau0
c     ----------------------------------------------
c
      wbas=igl**2.*3.
      wmin=exp(-(igl+0.5)**2/wbas)
c
      iml=nmdl
      itl=ntal
c
      dg1=abs(sazg(2)-sazg(1))
      dg2=abs(sozg(2)-sozg(1))
      dg3=abs(reag(2)-reag(1))
c
      rg0=(saz-sazg(1))/dg1+1.
      ig0=nint(rg0)
      do id0=-igl,igl
       ig1(id0)=ig0+id0
       if(ig1(id0).lt.1  ) ig1(id0)=1
       if(ig1(id0).gt.ng1) ig1(id0)=ng1
       rg1(id0)=ig1(id0)-rg0
      enddo
      rg0=(soz-sozg(1))/dg2+1.
      ig0=nint(rg0)
      do id0=-igl,igl
       ig2(id0)=ig0+id0
       if(ig2(id0).lt.1  ) ig2(id0)=1
       if(ig2(id0).gt.ng2) ig2(id0)=ng2
       rg2(id0)=ig2(id0)-rg0
      enddo
      rg0=(abs(rea)-reag(1))/dg3+1.
      ig0=nint(rg0)
      do id0=-igl,igl
       ig3(id0)=ig0+id0
       if(ig3(id0).lt.1  ) ig3(id0)=2
       if(ig3(id0).gt.ng3) ig3(id0)=ng3-1
       rg3(id0)=ig3(id0)-rg0
      enddo
c
      rr(1:nbsl,1:2)=0.
      ra(1:nbsl,1:itl,1:iml)=0.
      t0(1:nbsl,1:itl,1:iml)=0.
      t1(1:nbsl,1:itl,1:iml)=0.
      wgt=0.
      do id3=-igl,igl		! rea
       do id2=-igl,igl		! soz
        ix2=ix+ig2(id2)
        do id1=-igl,igl		! saz
         wt=exp(-(rg1(id1)**2+rg2(id2)**2+rg3(id3)**2)/wbas)-wmin
         if(wt.gt.0) then
          ix1=ix+ig1(id1)
          i=((ig3(id3)-1)*ng2+ig2(id2)-1)*ng1+ig1(id1)
          wgt=wgt+wt
          do nb=1,nbsl
           ra0=tbl(i,nb,1,1)
           rr(nb,1)=rr(nb,1)+tb0(i,nb)*wt
           rr(nb,2)=rr(nb,2)+ra0*wt
           do im=1,iml
            do it=1,itl
             t0(nb,it,im)=t0(nb,it,im)+tbl(ix2,nb,it,im)*wt
             t1(nb,it,im)=t1(nb,it,im)+tbl(ix1,nb,it,im)*wt
             ra(nb,it,im)=ra(nb,it,im)+(tbl(i,nb,it,im)-ra0)*wt
            enddo
           enddo
          enddo
         endif
        enddo
       enddo
      enddo
c
      do nb=1,nbsl
       rr(nb,1:2)=rr(nb,1:2)/wgt
       t0(nb,1:itl,1:iml)=t0(nb,1:itl,1:iml)/wgt
       t1(nb,1:itl,1:iml)=t1(nb,1:itl,1:iml)/wgt
       ra(nb,1:itl,1:iml)=ra(nb,1:itl,1:iml)/wgt
       ta(nb,1:itl,1:iml)=tbl(ngx+2,nb,1:itl,1:iml)	! AOT
       tb(nb,1:itl,1:iml)=tbl(ngx+1,nb,1:itl,1:iml)	! TOT
       sa(nb,1:itl,1:iml)=tbl(ngd,nb,1:itl,1:iml)
c      --------------------------- altitude correction
       rr(nb,2)=rr(nb,2)+rr(nb,1)*(epr1-1.)
       rr(nb,1)=rr(nb,1)*epr1
       rt0=t0(nb,1,1)**(epr1-1.)
       rt1=t1(nb,1,1)**(epr1-1.)
       tau0=tb(nb,1,1)
       t0(nb,1:itl,1:iml)=t0(nb,1:itl,1:iml)*rt0
       t1(nb,1:itl,1:iml)=t1(nb,1:itl,1:iml)*rt1
       sa(nb,1:itl,1:iml)=sa(nb,1:itl,1:iml)*epr1
       tb(nb,1:itl,1:iml)=tb(nb,1:itl,1:iml)-tau0+tau0*epr1
      enddo
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE read_ANC(iopr,iyr,jdy,nl,ml,nal, r4anc,ifanc,
     &                    rlon0,rlat0,rhr0,elv, tpw,atp,win,epr)
c
      implicit none
c
      integer	:: iopr,iyr,jdy,nl,ml,nal,ifanc
      real(4)	:: r4anc(nl*ml,nal,4)
      real(4)	:: rlon0,rlat0,rhr0,elv
      real(4)	:: tpw,atp,win,epr
c
      integer,parameter :: isk=23
      real(4)	:: r4buf(isk+nl*ml,2)
      character	:: fln*80
      character(len=10)	:: cdate
      character(len=19)	:: indirA
      character(len=19)	:: indirF
      integer	:: ifg1,ifg2
      real(4)	:: ranc(4)
      real(4)	:: wts(4),wtt(2)
      integer	:: nmi(4),ihr(2)
      integer	:: mnl(12)=(/0,31,59,90,120,151,181,212,243,273,304,334/)
      real(4)	:: def(4)=(/ 10.0, 288.15, 7.0, 101325./)
      integer	:: nml,nrec,iyl,j,iy1,jd,im,mnl1,im1,id1,ih1
      integer	:: np,nm,m1,m2,n1,n2,it,is
      real(4)	:: rdy,dd,rex,wta,avg,slp
      real(4)	:: rlon,rlat,rhr,rh,rm,rn
c
      nml=nl*ml
c
      if(iopr.le.0) then
c
       nrec=4*(isk+nml)
c
       if(mod(iyr,4).eq.0) then
        iyl=366
       else
        iyl=365
       endif
c
       do j=1,nal
        rdy=jdy+(j-1.)/(nal-1.)
        if(rdy.ge.iyl+1) then
         iy1=iyr+1
         rdy=rdy-iyl
        else
         iy1=iyr
        endif
c
        jd=int(rdy)
        do im=1,12
         mnl1=mnl(im)
         if((mod(iy1,4).eq.0).and.(im.gt.2)) mnl1=mnl1+1
         if(jd.gt.mnl1) then
          im1=im
          id1=jd-mnl1
         endif
        enddo
        ih1=nint((rdy-jd)*24)
        write(cdate,'(i4.4,3i2.2)') iy1,im1,id1,ih1
        write(indirA,'(a12,i4.4,i2.2,a1)') '../anc/GGLA/',iy1,im1,'/'
        write(indirF,'(a12,i4.4,i2.2,a1)') '../anc/GGLF/',iy1,im1,'/'
c
        do np=1,4
c
         ifg2=-1
         if(np.eq.1) then
          fln=indirA//'GGLA_'//cdate//'0000.PTW'
          open(14,file=fln,iostat=ifg1,
     &    access='direct',recl=nrec,form='unformatted',status='old')
          if(ifg1.ne.0) then
           fln=indirF//'GGLF_'//cdate//'0000.PTW'
           open(14,file=fln,iostat=ifg1,
     &     access='direct',recl=nrec,form='unformatted',status='old')
           ifanc=1
          else
           ifanc=0
          endif
c
         elseif(np.eq.2) then
          fln=indirA//'GGLA_'//cdate//'0000.T'
          open(14,file=fln,iostat=ifg1,
     &    access='direct',recl=nrec,form='unformatted',status='old')
          if(ifg1.ne.0) then
           fln=indirF//'GGLF_'//cdate//'0000.T'
           open(14,file=fln,iostat=ifg1,
     &     access='direct',recl=nrec,form='unformatted',status='old')
          endif
c
         elseif(np.eq.3) then
          fln=indirA//'GGLA_'//cdate//'0000.U'
          open(14,file=fln,iostat=ifg1,
     &    access='direct',recl=nrec,form='unformatted',status='old')
          if(ifg1.ne.0) then
           fln=indirF//'GGLF_'//cdate//'0000.U'
           open(14,file=fln,iostat=ifg1,
     &     access='direct',recl=nrec,form='unformatted',status='old')
          endif
          fln=indirA//'GGLA_'//cdate//'0000.V'
          open(15,file=fln,iostat=ifg2,
     &    access='direct',recl=nrec,form='unformatted',status='old')
          if(ifg2.ne.0) then
           fln=indirF//'GGLF_'//cdate//'0000.V'
           open(15,file=fln,iostat=ifg2,
     &     access='direct',recl=nrec,form='unformatted',status='old')
          endif
c
         elseif(np.eq.4) then
          fln=indirA//'GGLA_'//cdate//'0000.Psea'
          open(14,file=fln,iostat=ifg1,
     &    access='direct',recl=nrec,form='unformatted',status='old')
          if(ifg1.ne.0) then
           fln=indirF//'GGLF_'//cdate//'0000.Psea'
           open(14,file=fln,iostat=ifg1,
     &     access='direct',recl=nrec,form='unformatted',status='old')
          endif
         endif
c
         if(ifg1.eq.0) then
          if(np.eq.1) write(6,'(a)') 'Read '//trim(fln)//' for '//cdate
          read(14,rec=1,iostat=ifg1) (r4buf(nm,1),nm=1,isk+nml)
          close(14)
          if(ifg2.eq.0) then
           read(15,rec=1,iostat=ifg2) (r4buf(nm,2),nm=1,isk+nml)
           close(15)
          endif
         endif
c
         if(ifg2.eq.0) then
          do nm=isk+1,isk+nml
           r4buf(nm,1)=sqrt(r4buf(nm,1)**2+r4buf(nm,2)**2)
          enddo
         endif
c
         if(ifg1.ne.0) then
          if(np.eq.1) write(6,'(a)') 'No ANC data : '//trim(fln)
          do nm=1,nml
           r4anc(nm,j,np)=def(np)
          enddo
          ifanc=-1
         else
          do nm=1,nml
           r4anc(nm,j,np)=r4buf(isk+nm,1)
          enddo
         endif
c
        enddo	! np
c
       enddo	! nal
c
c     ---------------------------------------------------------- interp
      else
       dd=360./nl
       rex=5.25539	! (g=0.9807/(R=8.314/M=0.02896)/0.0065)
       rlon=rlon0
       rlat=rlat0
       rhr=rhr0
c
       rh=rhr/24.*(nal-1.)+1
       ihr(1)=int(rh)
       if(ihr(1).lt.1) then
        ihr(1)=1
       elseif(ihr(1).ge.nal) then
        ihr(1)=nal-1
       endif
       ihr(2)=ihr(1)+1
       wtt(1)=ihr(2)-rh
       wtt(2)=1.-wtt(1)
c
       rm=(90.-rlat)/dd+1.
       m1=int(rm)
       rm=rm-m1
       if(m1.lt. 1) m1=1
       if(m1.ge.ml) m1=ml-1
       m2=m1+1
       rn=(rlon-0.)/dd+1.
       n1=int(rn)
       rn=rn-n1
       n2=n1+1
       if(n1.lt. 1) n1=n1+nl
       if(n1.gt.nl) n1=n1-nl
       if(n2.lt. 1) n2=n2+nl
       if(n2.gt.nl) n2=n2-nl
       wts(1)=(1.-rn)*(1.-rm)
       wts(2)=(   rn)*(1.-rm)
       wts(3)=(1.-rn)*(   rm)
       wts(4)=(   rn)*(   rm)
       nmi(1)=(m1-1)*nl+n1
       nmi(2)=(m1-1)*nl+n2
       nmi(3)=(m2-1)*nl+n1
       nmi(4)=(m2-1)*nl+n2
c
       do np=1,4
        wta=0.
        avg=0.
        do it=1,2
         do is=1,4
          if((np.eq.3).or.(r4anc(nmi(is),ihr(it),np).ge.0)) then
           wta=wta+wts(is)*wtt(it)
           avg=avg+r4anc(nmi(is),ihr(it),np)*wts(is)*wtt(it)
          endif
         enddo
        enddo
        if(wta.gt.0) then
         ranc(np)=avg/wta
         if(ranc(np).lt.0.) ranc(np)=0.
        else
         ranc(np)=def(np)
        endif
       enddo
c
       tpw=ranc(1)
       atp=ranc(2)
       win=ranc(3)
       slp=ranc(4)
       if(elv.le.0.) then
        epr=slp/101325.
       else
        epr=slp/101325.*(atp/(atp+0.0065*elv))**rex
       endif
c
      endif
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE read_TOZ(iopr,iyr,jdy,ndyl, i2ozn,ifozn,
     &                    rlon,rlat,rhr, ozn)
c
      implicit none
c
      character(len=6) :: indir='../anc'
      real(4),parameter	:: toz0=343.79
      integer,parameter :: nl1=360,ml1=180
c
      integer	:: iopr,iyr,jdy,ndyl
      integer*2 :: i2ozn(nl1*ml1,ndyl)
      integer	:: ifozn
      real(4)	:: rlon,rlat,rhr
      real(4)	:: ozn
c
      character(len=200) :: oznf
      integer*2 :: i2buf(288,145)		! buf for TOVS
      real(4)	:: wtt(2),rdy
      integer	:: ihr(2)
      integer	:: ierr
      integer	:: ndy,jdy1,ibk,jd1,iy1,im1,id1,im,id0
      integer	:: nl,ml,isk,m,n,m1,n1,nm,j
      integer	:: iw,md,mm,nd,nn
      integer	:: i2org,i2big,i2lit
      integer	:: ifoz0
      real(4)	:: dd1,rm,rn
      real(4)	:: toz2,wta2,toz,wta,wbas,wmin,wt
      integer	:: mnl(12)=(/0,31,59,90,120,151,181,212,243,273,304,334/)
c
      if(iopr.eq.0) then
       ifozn=0
       iy1=iyr
c
       do ndy=1,ndyl
        jdy1=jdy-ndy+2
c
c       ------------------------------------- search back dates
        ibk=0
        do while(ibk.le.10)
         jd1=jdy1-ibk
         im1=0
         if(jd1.lt.1) then
          im1=12
          iy1=iyr-1
          id1=jd1+31
         else
          do im=1,12
           if((mod(iy1,4).eq.0).and.(im.gt.2)) then
            id0=jd1-mnl(im)-1
           else
            id0=jd1-mnl(im)
           endif
           if(id0.ge.1) then
            id1=id0
            im1=im
           endif
          enddo
         endif
         nl=nl1
         ml=ml1
         isk=1
c
         write(oznf,'(a11,i4.4,i2.2,a6,i4.4,2i2.2,a17)')
     &   '/ozone_jma/',iy1,im1,'/OZONE',iy1,im1,id1,'120000.ctm_le.ana'
         open(19,file=indir//oznf,iostat=ierr,
     &   access='direct',recl=2*288,form='unformatted',status='old')
         if(ierr.ne.0) then
          write(oznf,'(a7,i4.4,i2.2,a6,i4.4,2i2.2,a17)')
     &    '/ozone/',iy1,im1,'/OZONE',iy1,im1,id1,'120000.ctm_le.ana'
          open(19,file=indir//oznf,iostat=ierr,
     &    access='direct',recl=2*288,form='unformatted',status='old')
         endif
         if(ierr.eq.0) then
          nl=288
          ml=145
          isk=1
          ibk=999
          ifoz0=0
         endif
c
         if(ierr.ne.0) then
          write(oznf,'(a11,i4.4,i2.2,a6,i4.4,2i2.2,a17)')
     &    '/ozone_jma/',iy1,im1,'/OZONE',iy1,im1,id1,'120000.ctm_le.fcs'
          open(19,file=indir//oznf,iostat=ierr,
     &    access='direct',recl=2*288,form='unformatted',status='old')
          if(ierr.ne.0) then
           write(oznf,'(a7,i4.4,i2.2,a6,i4.4,2i2.2,a17)')
     &     '/ozone/',iy1,im1,'/OZONE',iy1,im1,id1,'120000.ctm_le.fcs'
           open(19,file=indir//oznf,iostat=ierr,
     &     access='direct',recl=2*288,form='unformatted',status='old')
          endif
          if(ierr.eq.0) then
           nl=288
           ml=145
           isk=1
           ibk=999
           ifoz0=1
          endif
         endif
c
         if(ierr.ne.0) then
          write(oznf,'(a7,i4.4,i2.2,a6,i4.4,2i2.2,a7)')
     &                '/ozone/',iy1,im1,'/OZONE',iy1,im1,id1,'.omi_le'
          open(19,file=indir//oznf,iostat=ierr,
     &    access='direct',recl=2*nl1,form='unformatted',status='old')
          if(ierr.eq.0) then
           ibk=999
           ifoz0=2
          endif
         endif
c
         if(ierr.ne.0) then
          write(oznf,'(a7,i4.4,i2.2,a6,i4.4,2i2.2,a7)')
     &          '/ozone/',iy1,im1,'/OZONE',iy1,im1,id1,'.epc_le'
          open(19,file=indir//oznf,iostat=ierr,
     &    access='direct',recl=2*nl1,form='unformatted',status='old')
          if(ierr.eq.0) then
           ibk=999
           ifoz0=2
          endif
         endif
c
         if(ierr.ne.0) then
          write(oznf,'(a7,i4.4,i2.2,a6,i4.4,2i2.2,a11)')
     &          '/ozone/',iy1,im1,'/OZONE',iy1,im1,id1,'.omi_le.nrt'
          open(19,file=indir//oznf,iostat=ierr,
     &    access='direct',recl=2*nl1,form='unformatted',status='old')
          if(ierr.eq.0) then
           ibk=999
           ifoz0=3
          endif
         endif
c
         if(ierr.ne.0) then
          write(oznf,'(a7,i4.4,i2.2,a8,i4.4,2i2.2)')
     &          '/ozone/',iy1,im1,'/TOVS_OZ',iy1,im1,id1
          open(19,file=indir//oznf,iostat=ierr,
     &    access='direct',recl=2*288,form='unformatted',status='old')
          if(ierr.eq.0) then
           nl=288
           ml=145
           isk=0
           ibk=999
           ifoz0=4
          endif
         endif
c
         ibk=ibk+1
        enddo
c       ----------------------------------------------
c
        if(ibk.eq.1000) then
         if(nl.eq.nl1) then
          do m=1,ml
           read(19,rec=isk+m,iostat=ierr) (i2ozn((m-1)*nl+n,ndy),n=1,nl)
          enddo
          close(19)
c
         elseif(ifoz0.ne.4) then
          do m1=1,ml
           read(19,rec=isk+m1,iostat=ierr) (i2buf(n1,m1),n1=1,288)
          enddo
          close(19)
          do m=1,ml1
           m1=nint((m-0.5)/1.25)+1
           do n=1,nl1
            n1=nint(144+(n-0.5)/1.25)+1
            if(n1.gt.288) n1=n1-288
            i2ozn((m-1)*nl1+n,ndy)=i2buf(n1,m1)
           enddo
          enddo
c
         else
          do m1=1,ml
           read(19,rec=isk+m1,iostat=ierr) (i2buf(n1,m1),n1=1,288)
          enddo
          close(19)
          do m=1,ml1
           m1=nint((m-0.5)/1.25)+1
           do n=1,nl1
            n1=nint((n-0.5)/1.25)+1
            i2org=i2buf(n1,m1)
            i2big=i2org/256
            i2lit=i2org-i2big*256
            i2ozn((m-1)*nl1+n,ndy)=(i2big+i2lit*256)*10
           enddo
          enddo
         endif
        else
         ierr=-1
         ifoz0=-1
        endif
c
        if(ierr.eq.0) then
         write(6,'(a)') 'Read '//indir//trim(oznf)
        else
         write(6,'(a33,i3)') ' No ozone data: = 343.79 DU, day:',ndy
         do nm=1,ml1*nl1
          i2ozn(nm,ndy)=nint(toz0*10)
         enddo
        endif
c
        if(ndy.eq.2) ifozn=ifoz0
       enddo
c
c     ------------------------------------------------------------------
      else
       dd1=360./nl1
c
       rdy=(1.0-rhr/24.)+1.5
       ihr(1)=int(rdy)
       if(ihr(1).lt.1) then
        ihr(1)=1
       elseif(ihr(1).ge.ndyl) then
        ihr(1)=ndyl-1
       endif
       ihr(2)=ihr(1)+1
       wtt(1)=ihr(2)-rdy
       wtt(2)=1.-wtt(1)
c
       rm=(89.50-rlat)/dd1+1.
       m1=nint(rm)
       rn=(179.50+rlon)/dd1+1.
       n1=nint(rn)
c
       toz2=0.
       wta2=0.
       do j=1,2
        ndy=ihr(j)
        iw=0
        do while (iw.lt.5)
         iw=iw+1
         toz=0.
         wta=0.
         wbas=iw**2.*3.
         wmin=exp(-(iw+0.5)**2/wbas)
         do md=-iw,iw
          mm=m1+md
          if((mm.gt.0).and.(mm.le.ml1)) then
           do nd=-iw,iw
            nn=n1+nd
            if(nn.lt.  1) nn=nn+nl1
            if(nn.gt.nl1) nn=nn-nl1
            wt=exp(-((nn-rn)**2+(mm-rm)**2)/wbas)-wmin
            nm=(mm-1)*nl1+nn
            if((i2ozn(nm,ndy).gt.0).and.(i2ozn(nm,ndy).lt.6500).and.
     &         (wt.gt.0)) then
             toz=toz+i2ozn(nm,ndy)*0.1*wt		! 0.1: slope
             wta=wta+wt
            endif
           enddo
          endif
         enddo
         if(wta.gt.0.1) then
          toz2=toz2+(toz/wta)*wtt(j)
          wta2=wta2+wtt(j)
          iw=999
         endif
        enddo
       enddo
       if(wta2.gt.0.) then
        ozn=toz2/wta2
       else
        ozn=toz0
       endif
c
      endif
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE SUNGLINT(TH0,TH,FI,WIND,nbul,nbu, sgr)
c
      implicit none
c
      real(4),parameter :: PI=3.14159265
      real(4),parameter :: RD=PI/180.
      real(4) :: cm(15)=
     &      (/1.3395,1.3383,1.3371,1.3351,1.3336,1.3327,1.3310,1.3310,
     &        1.3298,1.3287,1.3287,1.3267,1.3211,1.3168,1.2953/)
c
      real(4) :: TH0,TH,FI,WIND
      integer :: nbul,nbu(nbul)
      real(4) :: sgr(nbul)
c
      real(4) :: sth0,cth0,sth,cth,cdph,cos_w,s2,theta,sh,pw,y,z,buf
      integer :: nb
c
      if(TH0.lt.85.) then
       sth0=sin(TH0*RD)
       cth0=cos(TH0*RD)
       sth =sin(TH*RD)
       cth =cos(TH*RD)
       cdph=cos(FI*RD)
c
       buf=cth*cth0+sth*sth0*cdph
       if(buf.gt.1.) buf=1.
       cos_w=cos(acos(buf)/2.)
       buf=(cth+cth0)/(2.*cos_w)
       if(buf.gt.1.) buf=1.
       theta=acos(buf)
       sh=theta/RD
       if(sh.lt.85.) then
        s2=0.003+0.00512*WIND
        pw=exp(-tan(theta)**2./s2)/s2/(4.*cth0*cth*cos(theta)**4.)
        do nb=1,nbul
         buf=cm(nbu(nb))
         y=sqrt(buf**2.+cos_w**2.-1.0)/buf
         z=1./(cos_w+y*buf)**2.+1./(y+buf*cos_w)**2.
         sgr(nb)=pw*(1.-2.*buf*y*z*cos_w)
        enddo
       else
        sgr(1:nbul)=0.
       endif
c
      else
       sgr(1:nbul)=0.
      endif
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE estimate_aot(nmdl,nbl,nbu,wl,ra,t0,t1,sa,ta,tb,tg,
     &           rrc,rsc,sgr,amss,nbal,nba,ifngc, aot,rw)
c
      implicit none
c
      integer,parameter :: nbl0=15
      integer,parameter :: ntal=6
      integer,parameter :: nb1=7
      integer,parameter :: nb2=8
      integer,parameter :: nbtl=4
      integer	:: nbt(nbtl)=(/2,3,4,5/)	! in nbu
      real,parameter :: rcmin=0.0001	! rc min
      real,parameter :: rwmin=0.0002	! rw min at VIS
      real	:: pp(11)=(/-0.005,-0.363, 0.661,-1.078, 1.005,-1.426,
     &                  -2.112, 1.709, 2.407,-2.453, 2.066/)
c
c     -------------------- input
      integer :: nmdl,nbl
      integer :: nbu(nbl)
      real(4) :: wl(nbl)
      real(4) :: ra(nbl,ntal,nmdl)
      real(4) :: t0(nbl,ntal,nmdl)
      real(4) :: t1(nbl,ntal,nmdl)
      real(4) :: sa(nbl,ntal,nmdl)
      real(4) :: ta(nbl,ntal,nmdl)
      real(4) :: tb(nbl,ntal,nmdl)
      real(4) :: tg(nbl)
      real(4) :: rrc(nbl)
      real(4) :: rsc(nbl)
      real(4) :: sgr
      real(4) :: amss
      integer :: nbal
      integer :: nba(nbal)	! in nbu
      integer :: ifngc
c     -------------------- output
      real(4) :: aot(nbl)
      real(4) :: rw(nbl)
c     --------------------
c
      integer :: k,kl
      real(4) :: tb1(nbl)
      real(4) :: rc1(nbl)
      real(4) :: rpp,tr,tr0,drc,wpp
      real(4) :: rw1(nbal)
      real(4) :: wbd(nbal)
      real(4) :: wl1,wl2
      real(4) :: rsim(nbl,nmdl)
      integer :: itau(2,nmdl)
      real(4) :: wtau(2,nmdl)
      integer :: imdl(2)
      real(4) :: wmdl(2)
      real(4) :: rmm
      real(4) :: r0,rwb,rd
c
      integer :: nb0,nb,nbll
      integer :: it,nt
      integer :: im,nm
      integer :: im1,im2,it1,it2,im1r
      real(4) :: wt1,wt2,wt1r
      real(4) :: r1,r2
      real(4) :: wt,wta,w1
      real(4) :: ra2,tt2,sp2,ta2,bb
c
      nbll=1+nbtl+nbal
      aot(1:nbl)=-9.999
      rw(1:nbl)=-9.999
c
c     ---------------------------------- error-range -> weight
      wl2=wl(nba(1))
      do nb=1,nbal
       nb0=nba(nb)
       wl1=wl(nb0)
       rw1(nb)=rsc(nb0)
       if(rw1(nb).lt.rwmin) rw1(nb)=rwmin
       if(wl1.gt.wl2) then
        rw1(nb)=0.005	! error increase of far wl
       endif
      enddo
      wbd(1:nbal)=(rw1(1)/rw1(1:nbal))**2
c      if(wbd(nbal).gt.wbd(2)) write(6,'(7f9.4)') rw1(1:nbal)
c
      if(sgr.gt.0.0002) then
       kl=2
      else
       kl=1
      endif
c
      do k=1,kl
c
c      ---------------------------------- alignment interp
       if(sgr.gt.0.0002) then
        if(k.eq.1) then
         it=3
         im=nmdl-1
         tb1(1:nbl)=tb(1:nbl,it,im)
         wpp=(sgr-0.0002)/0.01
         if(wpp.gt.1.) wpp=1.
        else
         do nb=1,nbl
          tb1(nb)=0.
          do nm=1,2
           im=imdl(nm)
           do nt=1,2
            it=itau(nt,im)
            wt=wtau(nt,im)*wmdl(nm)
            wta=wta+wt
            tb1(nb)=tb1(nb)+tb(nb,it,im)*wt
           enddo
          enddo
          tb1(nb)=tb1(nb)/wta
         enddo
        endif
        tr0=exp(-tb1(nb1)*amss)*tg(nb1)
        do nb=1,nbl
         nb0=nbu(nb)
         if(nb0.le.11) then
          rpp=(pp(nb0)-pp(nb1))/(pp(nb2)-pp(nb1))
         else
          rpp=(0.0-pp(nb1))/(pp(nb2)-pp(nb1))
         endif
         tr=exp(-tb1(nb)*amss)*tg(nb)
         drc=(rrc(nb2)-rrc(nb1))*rpp
         rc1(nb)=rrc(nb)-drc*tr/tr0*wpp
         if(rc1(nb).lt.rcmin) rc1(nb)=rcmin
        enddo
       else
        rc1(1:nbl)=rrc(1:nbl)
       endif
c
c      ---------------------------------- tau
       do im=1,nmdl
        nb0=nba(1)
        it1=0
        it=1
        do while(it.lt.ntal)
         it2=it+1
         r1=ra(nb0,it ,im)+rsc(nb0)*t0(nb0,it ,im)*t1(nb0,it ,im)
     &                   /(1.-sa(nb0,it ,im)*rsc(nb0))
         r2=ra(nb0,it2,im)+rsc(nb0)*t0(nb0,it2,im)*t1(nb0,it2,im)
     &                   /(1.-sa(nb0,it2,im)*rsc(nb0))
         rd=r2-r1
         if(abs(rd).gt.1.0E-8) then
          w1=(r2-rc1(nb0))/rd
          if(((it.eq.1).and.(w1.ge.0.)).or.
     &       ((it.eq.ntal-1).and.(w1.le.1.)).or.
     &       ((w1.ge.0.).and.(w1.le.1.))) then
           it1=it
           wt1=w1
           it=ntal
          endif
         endif
         it=it+1
        enddo
        if(it1.eq.0) return		! cannot find tau
        if(wt1.lt.-2.) return	! too high tau
        if((it1.eq.1).and.(wt1.gt.1)) wt1=1.
        it2=it1+1
        wt2=1.-wt1
        itau(1,im)=it1
        itau(2,im)=it2
        wtau(1,im)=wt1
        wtau(2,im)=wt2
c
        rsim(1:nbl,im)=0.
        do nb=1,nbll
         if(nb.eq.1) then
          nb0=6	! green
          r0=rc1(nb0)
         elseif(nb.le.nbtl+1) then
          nb0=nbt(nb-1)
          r0=rsc(nb0)*0.1+rwmin*0.9
          if((nb.eq.nbtl+1).and.(r0.lt.rwb)) r0=rwb
         else
          nb0=nba(nb-1-nbtl)
          r0=rsc(nb0)
         endif
         ra2=ra(nb0,it1,im)*wt1+ra(nb0,it2,im)*wt2
         tt2=t0(nb0,it1,im)*t1(nb0,it1,im)*wt1
     &      +t0(nb0,it2,im)*t1(nb0,it2,im)*wt2
         sp2=sa(nb0,it1,im)*wt1+sa(nb0,it2,im)*wt2
         if(nb.eq.1) then
          bb=tt2+(r0-ra2)*sp2
          if(bb.lt.0.1) bb=0.1
          rwb=(r0-ra2)/bb*0.494	! min blue/green(chla=131.24)
          if(rwb.lt.rwmin) rwb=rwmin
         else
          bb=1.-sp2*r0
          if(bb.lt.0.1) bb=0.1
          rsim(nb0,im)=ra2+r0*tt2/bb
         endif
        enddo
c
       enddo
c
c      -------------------------------- model selection
       rmm=0.
       wta=0.
       do nb=2,nbal
        nb0=nba(nb)
        im1=0
        im=1
        do while(im.lt.nmdl)
         r1=rsim(nb0,im)
         r2=rsim(nb0,im+1)
         rd=r2-r1
         if(abs(rd).gt.1.0E-8) then
          w1=(r2-rc1(nb0))/rd
          if(((im.eq.1).and.(w1.ge.0.)).or.
     &       ((im.eq.nmdl-1).and.(w1.le.1.)).or.
     &       ((w1.ge.0.).and.(w1.le.1.))) then
           im1=im
           wt1=w1
           im=nmdl
          endif
         endif
         im=im+1
        enddo
        if(im1.gt.0) then
         if(wt1.gt. 5.) wt1= 5.
         if(wt1.lt.-4.) wt1=-4.
         rmm=rmm+(im1*wt1+(im1+1)*(1.-wt1))*wbd(nb)
         wta=wta+wbd(nb)
        endif
       enddo
       if(wta.eq.0.) return
c
       rmm=rmm/wta
       im1=int(rmm)
       if(im1.lt.   1) im1=1
       if(im1.ge.nmdl) im1=nmdl-1
       im2=im1+1
       wt1=im2-rmm
       wt2=1.-wt1
c
c      ---------------------- avoid negative rw
       if(ifngc.gt.0) then
        do nb=1,nbtl
         nb0=nbt(nb)
         ra2=rsim(nb0,im1)*wt1+rsim(nb0,im2)*wt2
         if(rc1(nb0).lt.ra2) then
          im1r=0
          im=1
          do while(im.lt.nmdl)
           r1=rsim(nb0,im)
           r2=rsim(nb0,im+1)
           rd=r2-r1
           if(abs(rd).gt.1.0E-8) then
            w1=(r2-rc1(nb0))/rd
            if(((im.eq.1).and.(w1.ge.0.)).or.
     &         ((im.eq.nmdl-1).and.(w1.le.1.)).or.
     &         ((w1.ge.0.).and.(w1.le.1.))) then
             im1r=im
             wt1r=w1
             im=nmdl
            endif
           endif
           im=im+1
          enddo
          if(im1r.gt.0) then
           im1=im1r
           wt1=wt1r
           if(wt1.gt. 5.) wt1= 5.
           if(wt1.lt.-4.) wt1=-4.
           im2=im1+1
           wt2=1.-wt1
          endif
         endif
        enddo
       endif
c
c      --------------------------------
       imdl(1)=im1
       imdl(2)=im2
       wmdl(1)=wt1
       wmdl(2)=wt2
      enddo		! kl
c
c     -------------------------------- rw calculation
      do nb=1,nbl
       wta=0.
       ta2=0.
       ra2=0.
       tt2=0.
       sp2=0.
       do nm=1,2
        im=imdl(nm)
        do nt=1,2
         it=itau(nt,im)
         wt=wtau(nt,im)*wmdl(nm)
         wta=wta+wt
         ta2=ta2+ta(nb,it,im)*wt
         ra2=ra2+ra(nb,it,im)*wt
         tt2=tt2+t0(nb,it,im)*t1(nb,it,im)*wt
         sp2=sp2+sa(nb,it,im)*wt
        enddo
       enddo
       aot(nb)=ta2/wta
       ra2=ra2/wta
       tt2=tt2/wta
       sp2=sp2/wta
       bb=tt2+(rc1(nb)-ra2)*sp2
       if(bb.lt.0.1) bb=0.1
       rw(nb)=(rc1(nb)-ra2)/bb
c      write(6,'(i4,99f8.4)') nb,ta2,ra2,tt2,sp2,rc1(nb),rw(nb)
      enddo
c
      if((rw(6).lt.rwmin).and.(rw(7).gt.rw(6))) then	! recover
       rw(6)=rw(7)
      endif
c      write(6,'(20f8.4)') rc1(1:10)
c      write(6,'(20f8.4)') rsc(1:10)
c      write(6,'(20f8.4)') rw(1:10)
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE rw_sim_lci(nbul,rcr, rss)
c
c     rsb:490, rsg:565, rsr:672, rsn:868
c
      implicit none
c
      integer	:: nbul
      real		:: rcr(nbul)
      real		:: rss(nbul)
c
      integer	:: nb
      real		:: pi=3.14159265
      real		:: edge1(2)=(/-0.005,0.005/)
      real		:: edge2(2)=(/ 0.020,0.030/)
      real		:: lcimin1=-0.02
      real		:: lcimin2=-0.03
      real		:: lcc1(3)=(/1.00000,-0.70230,-0.28826/)
      real		:: lcc2(3)=(/1.00000,-0.39203,-0.59614/)
      real		:: rsb,rsg,rsr,rsn
      real		:: lci1,lci2,r
c
      real :: rwmax=0.10
      real :: rwmin(11)=(/0.00147, 0.00202, 0.00253, 0.00348, 0.00513,
     &           0.00274, 0.00022, 0.00022, 0.00002, 0.00001, 0.00001/)
c
      real :: sc(6,11)=reshape( (/
     &        0.00104, 0.09241, 0.00000,  0.00202,-0.28673, 0.00000,
     &        0.00141, 0.12084, 0.00000,  0.00249,-0.34851, 0.00000,
     &        0.00186, 0.14629, 0.00000,  0.00233,-0.36140,20.97795,
     &        0.00269, 0.21252, 0.00000,  0.00313,-0.34247,25.96427,
     &        0.00266, 0.35279, 0.00000,  0.00454,-0.10035, 0.00000,
     &        0.00274, 0.27420, 6.85505,  0.01316, 0.25119, 0.00000,
     &        0.00082, 0.08215, 2.05365,  0.00263, 0.39518, 0.00000,
     &        0.00082, 0.08209, 2.05215,  0.00262, 0.39530, 0.00000,
     &        0.00016, 0.01646, 0.41158, -0.00007, 0.14582, 0.00000,
     &        0.00005, 0.00529, 0.13221, -0.00010, 0.04316, 0.00000,
     &        0.00005, 0.00527, 0.13169, -0.00010, 0.04294, 0.00000
     &         /),(/6,11/) )
c
      rsb=rcr(4)
      if(rsb.lt.rwmin(4)) rsb=rwmin(4)
      rsg=rcr(6)
      if(rsg.lt.rwmin(6)) rsg=rwmin(6)
      rsr=rcr(7)
      if(rsr.lt.rwmin(7)) rsr=rwmin(7)
      rsn=rcr(10)
      if(rsn.lt.rwmin(10)) rsn=rwmin(10)
c
      lci1=lcc1(1)*rsg+lcc1(2)*rsb+lcc1(3)*rsn
      if(lci1.lt.lcimin1) lci1=lcimin1
      lci2=lcc2(1)*rsr+lcc2(2)*rsb+lcc2(3)*rsn
      if(lci2.lt.lcimin2) lci2=lcimin2
c
      do nb=1,nbul
       if(nb.gt.11) then
        rss(nb)=rwmin(11)
       else
        if(nb.ge.6) then		! wl>540
         r=(edge2(2)-lci1)/(edge2(2)-edge2(1))
        else
         r=(lci2-edge1(1))/(edge1(2)-edge1(1))
        endif
        if(r.gt.1.0) then
         rss(nb)=(sc(1,nb)+sc(2,nb)*lci1+sc(3,nb)*lci1**2)*pi
        elseif(r.lt.0.0) then
         rss(nb)=(sc(4,nb)+sc(5,nb)*lci2+sc(6,nb)*lci2**2)*pi
        else
         rss(nb)=((sc(1,nb)+sc(2,nb)*lci1+sc(3,nb)*lci1**2)*r
     &           +(sc(4,nb)+sc(5,nb)*lci2+sc(6,nb)*lci2**2)*(1.0-r))*pi
        endif
c
        if(rss(nb).gt.rwmax) then
         rss(nb)=rwmax
        elseif(rss(nb).lt.rwmin(nb)) then
         rss(nb)=rwmin(nb)
        endif
       endif
      enddo
c
      return
      end
c=======================================================================
c
c=======================================================================
      subroutine brdf_oc(ifg,nbwl,nbul,wl,saz,soz,rea,chl,win,
     &                   rgo0,foQ0, foQ1)
c
      implicit none
c
      integer,parameter :: N_W=7
      integer,parameter :: N_S=6
      integer,parameter :: N_C=6
      integer,parameter :: N_N=17
      integer,parameter :: N_A=13
      integer,parameter :: N_K=5
c
      real(4) :: WVL(N_W)=(/412.5,442.5,490.0,510.0,560.0,620.0,660.0/)
      real(4) :: SUN(N_S)=(/0,15,30,45,60,75/)
      real(4) :: LCHL(N_C)=(/-1.523,-1.000,-0.523, 0.000, 0.477, 1.000/)
      real(4) :: NADIR(N_N)=(/1.078, 3.411, 6.289, 9.278,12.300,15.330,
     &                       18.370,21.410,24.450,27.500,30.540,33.590,
     &                       36.640,39.690,42.730,45.780,48.26818296/)
      real(4) :: AZIM(N_A)=(/0,15,30,45,60,75,90,105,120,135,150,165,
     &                       180/)
      real(4) :: win0(5)=(/0,4,6,8,16/)
c
      integer :: ifg
      integer :: nbwl,nbul
      real(4) :: wl(nbul)
      real(4) :: saz,soz,rea,chl,win
      real(4) :: rgo0(N_N,N_K)
      real(4) :: foQ0(N_A,N_N,N_C,N_S,nbwl)
      real(4) :: foQ1(nbwl)
c
      real(4) :: foQ(N_A,N_N,N_C,N_S,N_W)
      real(4) :: foQn(nbwl)
      real(4) :: rgoth(91)
      real(4) :: rgo1
      real(4) :: rgon
      real(4) :: km=1.33
      real(4) :: d2r=3.141593/180.
      real(4) :: lch,sazw,rea1
      real(4) :: w0(0:1)
      real(4) :: w1(0:1)
      real(4) :: w2(0:1)
      real(4) :: w3(0:1)
      real(4) :: w4(0:1)
      real(4) :: wt
      integer :: i0,i1,i2,i3,i4,ierr,nb
      integer :: NW,NS,NC,NN,NA,NK
c
c     --------------------------------------------------------
      if(ifg.eq.1) then
       open(10,file='../lut/rgoth_Ebuchi.txt',iostat=ierr,
     & access='sequential',form='formatted',status='old')
c
       if(ierr.eq.0) then
        write(6,*) 'read = ../lut/rgoth_Ebuchi.txt'
        do i4=1,N_K
         do i1=1,91
          read(10,*) NK,NN,rgoth(i1)
         enddo
         do NN=1,N_N
          sazw=NADIR(NN)
          i1=1+int(sazw)
          w1(0)=i1-sazw
          w1(1)=1.-w1(0)
          rgo0(NN,i4)=rgoth(i1)*w1(0)+rgoth(i1+1)*w1(1)
         enddo
        enddo
        close(10)
       else
        rgo0(1:N_N,1:N_K)=0.528872
        write(6,*) 'No file: lut/rgoth_Ebuchi.txt'
       endif
c
c     ---------------
c       open(10,file='../lut/foq.dat',iostat=ierr,
       open(10,file='../lut/values0.dat',iostat=ierr,
     & access='sequential',form='formatted',status='old')
c
       if(ierr.eq.0) then
c        write(6,*) 'read = ../lut/foq.dat'
        write(6,*) 'read = ../lut/values0.dat'
        do NW=1,N_W
         do NS=1,N_S
          do NC=1,N_C
           read(10,*) wt		! for values0.dat
           do NN=1,N_N
c            read(10,'(f6.4,12f7.4)') (foQ(NA,NN,NC,NS,NW),NA=1,N_A)
            read(10,'(13f8.5)') (foQ(NA,NN,NC,NS,NW),NA=1,N_A)
           enddo
          enddo
         enddo
        enddo
        close(10)
       else
c        write(6,*) 'Cannot read: ../lut/foq.dat'
        write(6,*) 'Cannot read: ../lut/values0.dat'
        foQ(1:N_A,1:N_N,1:N_C,1:N_S,1:N_W)=1.0
       endif
c
       do nb=1,nbwl
        if(wl(nb).lt.WVL(1)) then
         foQ0(1:N_A,1:N_N,1:N_C,1:N_S,nb)=
     &                           foQ(1:N_A,1:N_N,1:N_C,1:N_S,1)
        elseif(wl(nb).gt.WVL(N_W)) then
         foQ0(1:N_A,1:N_N,1:N_C,1:N_S,nb)=
     &                           foQ(1:N_A,1:N_N,1:N_C,1:N_S,N_W)
        else
         do NW=1,N_W-1
          if((wl(nb).ge.WVL(NW)).and.(wl(nb).lt.WVL(NW+1))) then
           foQ0(1:N_A,1:N_N,1:N_C,1:N_S,nb)=
     &         (foQ(1:N_A,1:N_N,1:N_C,1:N_S,NW)*(WVL(NW+1)-wl(nb))
     &         +foQ(1:N_A,1:N_N,1:N_C,1:N_S,NW+1)*(wl(nb)-WVL(NW)))
     &         /(WVL(NW+1)-WVL(NW))
          endif
         enddo
        endif
       enddo
       ifg=2
      endif
c     --------------------------------------------------------
c
      lch=log10(chl)
      rea1=abs(180.-rea)
      sazw=asin(sin(saz*d2r)/km)/d2r
c
      if(soz.le.SUN(1)) then
       w0(0)=1.
       i0=1
      elseif(soz.ge.SUN(N_S)) then
       w0(0)=0.
       i0=N_S-1
      else
       do NS=1,N_S-1
        if((soz.ge.SUN(NS)).and.(soz.lt.SUN(NS+1))) then
         w0(0)=(SUN(NS+1)-soz)/(SUN(NS+1)-SUN(NS))
         i0=NS
        endif
       enddo
      endif
      w0(1)=1.-w0(0)
c
      if(sazw.le.NADIR(1)) then
       w1(0)=1.
       i1=1
      elseif(sazw.ge.NADIR(N_N)) then
       w1(0)=0.
       i1=N_N-1
      else
       do NN=1,N_N-1
        if((sazw.ge.NADIR(NN)).and.(sazw.lt.NADIR(NN+1))) then
         w1(0)=(NADIR(NN+1)-sazw)/(NADIR(NN+1)-NADIR(NN))
         i1=NN
        endif
       enddo
      endif
      w1(1)=1.-w1(0)
c
      if(rea1.le.AZIM(1)) then
       w2(0)=1.
       i2=1
      elseif(rea1.ge.AZIM(N_A)) then
       w2(0)=0.
       i2=N_A-1
      else
       do NA=1,N_A-1
        if((rea1.ge.AZIM(NA)).and.(rea1.lt.AZIM(NA+1))) then
         w2(0)=(AZIM(NA+1)-rea1)/(AZIM(NA+1)-AZIM(NA))
         i2=NA
        endif
       enddo
      endif
      w2(1)=1.-w2(0)
c
      if(lch.le.LCHL(1)) then
       w3(0)=1.
       i3=1
      elseif(lch.ge.LCHL(N_C)) then
       w3(0)=0.
       i3=N_C-1
      else
       do NC=1,N_C-1
        if((lch.ge.LCHL(NC)).and.(lch.lt.LCHL(NC+1))) then
         w3(0)=(LCHL(NC+1)-lch)/(LCHL(NC+1)-LCHL(NC))
         i3=NC
        endif
       enddo
      endif
      w3(1)=1.-w3(0)
c
      foQ1(1:nbwl)=0.
      foQn(1:nbwl)=0.
      do NS=0,1
       do NN=0,1
        do NA=0,1
         do NC=0,1
          wt=w0(NS)*w1(NN)*w2(NA)*w3(NC)
          do nb=1,nbwl
           foQ1(nb)=foQ1(nb)+foQ0(i2+NA,i1+NN,i3+NC,i0+NS,nb)*wt
           foQn(nb)=foQn(nb)+foQ0(1,1,i3+NC,1,nb)*wt
          enddo
         enddo
        enddo
       enddo
      enddo
      do nb=1,nbwl
       foQ1(nb)=foQ1(nb)/foQn(nb)
      enddo
c
c     -----
      if(win.le.win0(1)) then
       w4(0)=1.
       i4=1
      elseif(win.ge.win0(N_K)) then
       w4(0)=0.
       i4=N_K-1
      else
       do NK=1,N_K-1
        if((win.ge.win0(NK)).and.(win.lt.win0(NK+1))) then
         w4(0)=(win0(NK+1)-win)/(win0(NK+1)-win0(NK))
         i4=NK
        endif
       enddo
      endif
      w4(1)=1.-w4(0)
c
      rgo1=0.
      rgon=0.
      do NN=0,1
       do NK=0,1
        wt=w1(NN)*w4(NK)
        rgo1=rgo1+rgo0(i1+NN,i4+NK)*wt
        rgon=rgon+rgo0(i1+NN,1)*wt
       enddo
      enddo
      rgo1=rgo1/rgon
c
      foQ1(1:nbwl)=foQ1(1:nbwl)*rgo1
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE cal_iop(ifiop,rrs, chla,aph442,adg442,bbp442)
c
      implicit none
c
      integer,parameter :: ifocr=0	! 1:OC4 only, 0*:mix
      integer,parameter :: nbl0=7	! band number of rrs
      integer,parameter :: nbl=4	! bands for LMI
      real,parameter	:: rrsmin=0.0001
c
      integer :: ifiop
      REAL(4) :: rrs(nbl0)
      REAL(4) :: chla,aph442,adg442,bbp442
c
      REAL(4) :: nbc(nbl) =(/3,4,5,6/)
      REAL(4) :: rrsw(nbl)
c
      REAL(4) :: wvc(nbl0)=(/380.03, 412.51, 443.24, 489.85,
     &                       529.64, 566.16, 672.00/)
      REAL(4) :: aw(nbl0) =(/0.00377,0.00312,0.00510,0.01337, 
     &                       0.04211,0.06759,0.44573/)
      REAL(4) :: bbw(nbl0)=(/0.00472,0.00333,0.00239,0.00157, 
     &                       0.00112,0.00086,0.00041/)
      REAL(4) :: aph(nbl0)=(/0.65519,0.84176,0.98787,0.61767, 
     &                       0.30412,0.14724,0.58191/)
      REAL(4) :: adg(nbl0)=(/2.48084,1.54332,0.98481,0.49842, 
     &                       0.27961,0.16400,0.03501/)
      REAL(4) :: bbp(nbl0)=(/1.18327,1.08016,0.99723,0.89237, 
     &                       0.81833,0.75988,0.62823/)
      REAL(4) :: apg(nbl0)
c
      REAL(4) :: oc4(5)=(/0.39747,-3.42876, 5.33109,-5.39966, 1.73379/)
      REAL(4) :: oci(2)=(/-0.38817,236.59825/)
c
      integer :: nb,nb0,k,iw
      REAL(8) :: sxx,sxy,sxz,syy,syz,szz,x,y,z,ss,u
      REAL(8) :: v(nbl)
      REAL(4) :: rmax,r,ci,cx,wci,chli
      REAL(4) :: g1=0.0949
      REAL(4) :: g2=0.0794
      REAL(4) :: apg442,r1,r0,w,apg1,bbp1,df1,df0
c
      apg(1:nbl0)=aph(1:nbl0)*0.375+adg(1:nbl0)*(1.-0.375)
c
      chla  =-9.999
      apg442=-9.999
      bbp442=-9.999
      aph442=-9.999
      adg442=-9.999
c
c     ---------------- OCx & CI
      if(rrs(6).gt.rrsmin) then
       rmax=-9.999
       do nb=3,5
        if(rrs(nb).gt.rmax) rmax=rrs(nb)
       enddo
       if(rmax.gt.rrsmin) then
        r=log10(rmax/rrs(6))
        if(r.gt.1.55) r=1.55
c       if(r.gt.-0.3062) then	! chla<131.0539
        if(r.gt.-0.3275) then	! chla<200.3539
         chla=10.**(oc4(1)+oc4(2)*r+oc4(3)*r**2+oc4(4)*r**3+oc4(5)*r**4)
        endif
       endif
      endif
c
      if((ifocr.eq.0).and.(chla.gt.0.)) then
       if((rrs(3).gt.0.).and.(rrs(6).gt.0.).and.(rrs(7).gt.0.)) then
        ci=rrs(6)-(rrs(3)*(wvc(7)-wvc(6))+rrs(7)*(wvc(6)-wvc(3)))
     &           /(wvc(7)-wvc(3))
        cx=oci(1)+oci(2)*ci
        if(cx.lt.2.) then
         chli=10.**cx
         if(chli.lt.0.005) chli=0.005
         wci=((-0.0002)-ci)/((-0.0002)-(-0.0006))
         if(chla.lt.0.01) then
          wci=1.
         elseif(wci.gt.1) then
          wci=1.
         elseif(wci.lt.0) then
          wci=0.
         endif
         chla=chli*wci+chla*(1.-wci)
        endif
       endif
c      write(6,'(10f9.5)') chla,chli,ci,(rrs(nb),nb=1,7)
      endif
c
      if(ifiop.eq.1) then
c      ---------------- IOP by Linear Matrix Inversion 
       do nb=1,nbl
        nb0=nbc(nb)
        rrsw(nb)=rrs(nb0)/(0.529+1.70*rrs(nb0))	! under water
       enddo
c
       if((rrsw(2).gt.0.).and.(rrsw(3).gt.0.).and.(rrsw(4).gt.0.)) then
        sxx=0.
        sxy=0.
        sxz=0.
        syy=0.
        syz=0.
        szz=0.
        do nb=2,nbl
         nb0=nbc(nb)
         v(nb)=1-2*g2/(-g1+sqrt(g1**2+4*g2*rrsw(nb)))
         x=apg(nb0)
         y=bbp(nb0)*v(nb)
         z=-(aw(nb0)+bbw(nb0)*v(nb))
         sxx=sxx+x*x
         sxy=sxy+x*y
         sxz=sxz+x*z
         syy=syy+y*y
         syz=syz+y*z
         szz=szz+z*z
        enddo
        ss=sxx*syy-sxy*sxy
        if(abs(ss).gt.1D-10) then
         apg442=(sxz*syy-syz*sxy)/ss
         bbp442=(syz*sxx-sxz*sxy)/ss
c
         if(bbp442.lt.0.0001) then
          bbp1=0.0001
          apg442=0.
          do nb=2,nbl
           nb0=nbc(nb)
           apg442=apg442
     &           -((bbp(nb0)*bbp1+bbw(nb0))*v(nb)+aw(nb0))/apg(nb0)
          enddo
          apg442=apg442/(nbl-1)
          bbp442=-9.999
         else
          bbp1=bbp442
         endif
c
         if((apg442.gt.0.0001).and.(rrsw(1).gt.0.)) then
          df0=9999.
          r0=0.5
          do k=1,2
           do iw=-4,4
            if(k.eq.1) then
             r=r0+0.10*iw	! 0.1 ~ 0.9
            else
             r=r0+0.02*iw	! -0.08 ~ +0.08
            endif
            w=1.0+0.05*((r-0.5)/0.3)**4
            df1=0.
            do nb=1,3
             nb0=nbc(nb)
             apg1=(aph(nb0)*r+adg(nb0)*(1.-r))*apg442
             u=(bbw(nb0)+bbp(nb0)*bbp1)
     &        /(bbw(nb0)+bbp(nb0)*bbp1+aw(nb0)+apg1)
             df1=df1+(g1*u+g2*u**2-rrsw(nb))**2.*w
            enddo
            if(df1.lt.df0) then
             r1=r
             df0=df1
            endif
           enddo
           r0=r1
          enddo
          aph442=apg442*r0
          adg442=apg442*(1.-r0)
          bbp442=bbp1
         else
          apg442=-9.999
         endif
c
        endif
       endif
      endif
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE cal_par(nbl,nbpl,nbpi,iyr,jdy,rhr,rlat,csoz,rc,rs,
     &                   tgn,t0,t1,sa,f0w, par)
c
      implicit none
c
      REAL(8),parameter :: pi=3.14159265359
      REAL(8),parameter :: d2r=pi/180.
      integer,parameter :: nbv=4
      integer,parameter :: ifday=1
      real(4),parameter :: ah=10.	! <85deg to avoid too long amass
c
      integer :: nbl,nbpl
      integer :: nbpi(nbpl)
      integer :: iyr,jdy
      real(4) :: rhr,rlat,csoz,rc(nbl),rs(nbl)
      real(4) :: tgn(nbl),t0(nbl),t1(nbl),sa(nbl)
      real(4) :: f0w(nbpl)
      real(4) :: par
c
      real(4) :: rsoz,csozd1
      real(4) :: csozd(12),csoz1,trn(nbpl),trn0,trn1
      integer :: iyl(4)=(/0,366,731,1096/)
      real(4) :: dra,dse,rj,rjh,sig,sodec,omega
      real(4) :: wta,wtb,Acld,rs1,par1,par0,amass0
      integer :: ih,nb,nb0,nbc
c     ------------------------------------------------------
c
      dra=(0.9856003*jdy+357.528)*d2r
      dse=1.00014-0.01671*cos(dra)-0.00014*cos(2.0*dra)
c
      if(ifday.gt.0) then
       rj=int((iyr-2000)/4)*1461+iyl(mod(iyr,4)+1)+jdy
       rjh=rj+rhr/24.-0.5
       sig=rjh/365.25*2*pi
       sodec=0.006918-0.399912*cos(sig)+0.070257*sin(sig)
     &      -0.006758*cos(2*sig)+0.000907*sin(2*sig)
     &      -0.002697*cos(3*sig)+0.001480*sin(3*sig)
c
       do ih=1,12		! LT 0:30~11:30hr
        omega=(11.5+ih)/12.*pi
        csozd1=sin(sodec)*sin(rlat*d2r)
     &        +cos(sodec)*cos(rlat*d2r)*cos(omega)
        if(csozd1.lt.0.342) then
         rsoz=acos(csozd1)
         csozd(ih)=sin(rsoz)/(6378.+ah)*ah
     &            /sin(rsoz-asin(6378./(6378+ah)*sin(rsoz)))
        else
         csozd(ih)=csozd1
        endif
       enddo
      endif
c
      if(csoz.lt.0.342) then
       rsoz=acos(csoz)
       csoz1=sin(rsoz)/(6378.+ah)*ah
     &      /sin(rsoz-asin(6378./(6378+ah)*sin(rsoz)))
      else
       csoz1=csoz
      endif
c
      trn(1:nbpl)=-9.999
      trn0=0.
      wta=0.
      nbc=0
      do nb=1,nbpl
       nb0=nbpi(nb)
       if((rc(nb0).ge.0).and.(rs(nb0).ge.0).and.(nbc.lt.nbv)) then
        Acld=rc(nb0)/(t0(nb0)*t1(nb0)+rc(nb0)*sa(nb0))
        rs1=rs(nb0)
        if(rs1.gt.0.99) rs1=0.99
        if(Acld.lt.rs1) Acld=rs1
        if(Acld.gt.0.99) Acld=0.99
        trn(nb)=1.0/(1.-sa(nb0)*Acld)*(1.-Acld)/(1.-rs1)
        trn0=trn0+trn(nb)*f0w(nb)
        wta=wta+f0w(nb)
        nbc=nbc+1
       endif
      enddo
c
      if(wta.gt.0) then
       trn0=trn0/wta
       do nb=1,nbv
        if(trn(nb).lt.0.) trn(nb)=trn0
       enddo
c
       par1=0.
       wtb=0.
       do nb=1,nbv
        nb0=nbpi(nb)
        trn1=tgn(nb0)**(1./csoz1)*t0(nb0)*trn(nb)
        if(ifday.gt.0) then
         par0=0.
         do ih=1,12
          if(csozd(ih).gt.0.0175) then
           amass0=csoz1/csozd(ih)
           par0=par0+(trn1**amass0)*csozd(ih)
          endif
         enddo
         par0=par0/12.
        else
         par0=trn1*csoz
        endif
        par1=par1+par0*f0w(nb)
        wtb=wtb+f0w(nb)
       enddo
       if(ifday.gt.0) then
        par=(par1/wtb)*1743.21/dse**2 *0.1192	! Ein(mol)/m^2/day
       else
        par=(par1/wtb)*1743.21/dse**2 *0.1192/0.086400	! umol/m^2/s
       endif
c
      else
       par=-9.999
      endif
c
      return
      end
c=======================================================================
c
c=======================================================================
c     Read 2D int or real data to INT32 or FLOAT32
      SUBROUTINE read2d_full(file_id,dsetname,ml,nl, 
     &                       npix,nlin,i4buf,r4buf,slope,offst,error)
c
c     input: file_id,dsetname,ml,nl=>dim size
c     output: npix,nlin,(dset_id,dspace_id,mspace_id,class),i4buf or r4buf
c
      USE HDF5	! This module contains all necessary modules
c
      implicit none
c
      INTEGER(HID_T)	:: file_id		! File identifier
      CHARACTER(LEN=80)	:: dsetname		! Dataset name
      INTEGER			:: ml,nl,npix,nlin
c
      INTEGER(HID_T)	:: dset_id		! Dataset identifier
      INTEGER(HID_T)	:: type_id		! Datatype identifier
      INTEGER(HID_T)	:: dspace_id	! Dataspace identifier
      INTEGER(HID_T)	:: mspace_id	! memory identifier
      INTEGER(HID_T)	:: attr_id
      INTEGER(HSIZE_T)	:: adims(1)
      INTEGER(HSIZE_T)	:: dims(2)
      INTEGER(HSIZE_T)	:: maxdims(2)
      INTEGER(HSIZE_T)	:: start(2)
      INTEGER(HSIZE_T)	:: edge(2)
      CHARACTER(LEN=80)	:: attr_name	! Attribute name
      INTEGER			:: error0		! Error flag
      INTEGER			:: error		! Error flag
      INTEGER			:: rank
      INTEGER			:: class
      real(4)			:: slope,offst
c
      integer	::	i4buf(nl*ml)
      real(4)	::	r4buf(nl*ml)
c
      CALL h5eset_auto_f(0,error0)
      CALL h5dopen_f(file_id, dsetname, dset_id, error)
      if(error.lt.0) return
c
      attr_name='Slope_reflectance'
      CALL h5aopen_by_name_f(dset_id,'.', attr_name, attr_id, error)
      if(error.lt.0) then
       attr_name='Slope'
       CALL h5aopen_by_name_f(dset_id,'.', attr_name, attr_id, error)
      endif
      if(error.gt.-1) then
       CALL h5aget_type_f(attr_id, type_id, error)
       CALL h5aread_f(attr_id, type_id, slope, adims, error)
       CALL h5aclose_f(attr_id, error)
      else
       slope=1.
      endif
c
      attr_name='Offset_reflectance'
      CALL h5aopen_by_name_f(dset_id,'.', attr_name, attr_id, error)
      if(error.lt.0) then
       attr_name='Offset'
       CALL h5aopen_by_name_f(dset_id,'.', attr_name, attr_id, error)
      endif
      if(error.gt.-1) then
       CALL h5aget_type_f(attr_id, type_id, error)
       CALL h5aread_f(attr_id, type_id, offst, adims, error)
       CALL h5aclose_f(attr_id, error)
      else
       offst=0.
      endif
c
      CALL h5eset_auto_f(1,error0)
c
      CALL h5dget_space_f(dset_id, dspace_id, error)
      CALL h5sget_simple_extent_ndims_f(dspace_id, rank,error)
      if(rank.ne.2) then
       write(6,*) '[Error] '//trim(dsetname)//' rank = ',rank
       call exit(7)
       stop
      endif
c
      CALL h5dget_type_f(dset_id, type_id,error)
      CALL h5tget_class_f(type_id, class,error)
      CALL h5sget_simple_extent_dims_f(dspace_id, dims,maxdims,error)
      npix=dims(1)
      nlin=dims(2)
c
      if((npix.le.nl).and.(nlin.le.ml)) then
       if    (class.eq.H5T_INTEGER_F) then
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, i4buf, dims,error)
       elseif(class.eq.H5T_FLOAT_F) then
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, r4buf, dims,error)
       else
        write(6,*) '[Error] '//trim(dsetname)//' class = ',class
        call exit(7)
        stop
       endif
      else
c      ----------------------------------- Avoid error
       if(npix.gt.nl) then
        write(6,'(a15,i7,a3,i7)') '[warning] npix:',npix,' > ',nl
        write(6,*) trim(dsetname) 
        dims(1)=nl
        npix=nl
       endif
       if(nlin.gt.ml) then
        write(6,'(a15,i7,a3,i7)') '[warning] nlin:',nlin,' > ',ml
        write(6,*) trim(dsetname) 
        dims(2)=ml
        nlin=ml
       endif
c
       start(1)=0
       start(2)=0
       edge(1) =npix
       edge(2) =nlin
       CALL h5screate_simple_f(rank,edge, mspace_id,error)
       CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F,
     &                            start, edge, error)
       if    (class.eq.H5T_INTEGER_F) then
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, i4buf, edge,
     &                 error,mspace_id,dspace_id)
       elseif(class.eq.H5T_FLOAT_F) then
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, r4buf, edge,
     &                 error,mspace_id,dspace_id)
       endif
       CALL h5sclose_f(mspace_id, error)
c      -----------------------------------
      endif
c
      CALL h5sclose_f(dspace_id, error)
      CALL h5dclose_f(dset_id, error)
c
      return
      end subroutine read2d_full
c=======================================================================
c
c=======================================================================
c     Read double data to FLOAT64
      SUBROUTINE h5read_r8(file_id,dsetname,nml, r8buf,error)
c
      USE HDF5	! This module contains all necessary modules
c
      INTEGER(HID_T)	:: file_id		! File identifier
      CHARACTER(LEN=80)	:: dsetname		! Dataset name
      INTEGER			:: nml
      real(8)			:: r8buf(nml)
      INTEGER			:: error		! Error flag
c
      INTEGER(HID_T)	:: dset_id		! Dataset identifier
      INTEGER(HID_T)	:: dspace_id	! Dataspace identifier
c      INTEGER(HID_T)	:: mspace_id	! memory identifier
      INTEGER(HSIZE_T)	:: dims(6)
      INTEGER(HSIZE_T)	:: maxdims(6)
c      INTEGER(HSIZE_T)	:: start(6)
c      INTEGER(HSIZE_T)	:: edge(6)
      INTEGER			:: error0		! Error flag
      INTEGER			:: rank
      INTEGER			:: i,nml0
c
      CALL h5eset_auto_f(0,error0)
c
      CALL h5dopen_f(file_id, dsetname, dset_id, error)
      if(error.lt.0) then
       write(6,*) 'No data (r8) : '//trim(dsetname)
       return
      endif
c
      CALL h5dget_space_f(dset_id, dspace_id, error)
      CALL h5sget_simple_extent_ndims_f(dspace_id, rank,error)
      CALL h5sget_simple_extent_dims_f(dspace_id, dims,maxdims,error)
      nml0=1
      do i=1,rank
       nml0=nml0*dims(i)
      enddo
c      write(6,*) trim(dsetname)
c      write(6,'(a7,8i8)') 'dims = ',(dims(i),i=1,rank)
      if(nml0.gt.nml) then
       write(6,*) '[Error] '//trim(dsetname)
       write(6,'(a7,i8,a3,i8)') 'r8 dim:',nml0,' > ',nml
       write(6,'(a7,6i8)') 'dims = ',(dims(i),i=1,rank)
       call exit(8)
      endif
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, r8buf, dims,error)
c
      CALL h5eset_auto_f(1,error0)
      CALL h5sclose_f(dspace_id, error)
      CALL h5dclose_f(dset_id, error)
c
      return
      end subroutine h5read_r8
c=======================================================================
c
c=======================================================================
c     Read data attribute
      SUBROUTINE read_dattr(dset_id,attr_name,ifg,
     &                      class,c1buf,i4para,r4para,error)
c
      USE HDF5	! This module contains all necessary modules
c
      implicit none
c
      INTEGER(HID_T)	:: dset_id
      CHARACTER(LEN=80)	:: attr_name	! Attribute name
      INTEGER			:: ifg
      INTEGER			:: class
      character(len=400):: c1buf
      integer			:: i4para
      real(4)			:: r4para
      INTEGER 			:: error 		! Error flag
c
      INTEGER(HID_T)	:: attr_id
      INTEGER(HID_T)	:: type_id
      INTEGER(SIZE_T)	:: precision
      INTEGER(HSIZE_T)	:: dims(5)
      LOGICAL			:: f_corder_valid
      INTEGER			:: corder,cset
      INTEGER(HSIZE_T)	:: attrlen		! Length of the attribute string
      INTEGER(SIZE_T)	:: anmlen=50
      INTEGER			:: attr_num
      INTEGER(HSIZE_T)	:: attr_idx
      INTEGER			:: ierr
      real*8			:: r8para
c
      CALL h5eset_auto_f(0,ierr)
c
      if(attr_name.ne.'') then
       CALL h5aopen_by_name_f(dset_id,'.', attr_name, attr_id, error)
       if(error.gt.-1) then
        CALL h5aget_info_f(attr_id,f_corder_valid,corder,cset,attrlen,
     &                     error)
        CALL h5aget_type_f(attr_id, type_id, error)
        CALL h5tget_class_f(type_id, class,error)
        CALL h5tget_precision_f(type_id, precision  ,error)
        if(class.eq.H5T_INTEGER_F) then
         CALL h5aread_f(attr_id, type_id, i4para, dims, error)
         if(ifg.gt.0) write(6,*) attr_name,'=',i4para
        elseif(class.eq.H5T_FLOAT_F) then
         if(precision.eq.32) then
          CALL h5aread_f(attr_id, type_id, r4para, dims, error)
         else
          CALL h5aread_f(attr_id, type_id, r8para, dims, error)
          r4para=r8para
         endif
         if(ifg.gt.0) write(6,*) attr_name,'=',r4para
        elseif(class.eq.H5T_STRING_F) then
         CALL h5aread_f(attr_id, type_id, c1buf, dims, error)
         if(ifg.gt.0) write(6,*) attr_name,'=',c1buf(1:attrlen)
        endif
        CALL h5aclose_f(attr_id, error)
       endif
c
      else
       CALL h5aget_num_attrs_f(dset_id, attr_num, error)
       do attr_idx=0,attr_num-1
        CALL h5aopen_by_idx_f(dset_id, '.', H5_INDEX_CRT_ORDER_F,
     &       H5_ITER_INC_F, attr_idx, attr_id, error)
        CALL h5aget_name_f(attr_id, anmlen, attr_name, error)
        CALL h5aget_info_f(attr_id,f_corder_valid,corder,cset,attrlen,
     &                     error)
        CALL h5aget_type_f(attr_id, type_id, error)
        CALL h5tget_class_f(type_id, class,error)
        CALL h5tget_precision_f(type_id, precision  ,error)
        if(class.eq.H5T_INTEGER_F) then
         CALL h5aread_f(attr_id, type_id, i4para, dims, error)
         if(ifg.gt.0) write(6,*) attr_name,'=',i4para
        elseif(class.eq.H5T_FLOAT_F) then
         if(precision.eq.32) then
          CALL h5aread_f(attr_id, type_id, r4para, dims, error)
         else
          CALL h5aread_f(attr_id, type_id, r8para, dims, error)
          r4para=r8para
         endif
         if(ifg.gt.0) write(6,*) attr_name,'=',r4para
        elseif(class.eq.H5T_STRING_F) then
         CALL h5aread_f(attr_id, type_id, c1buf, dims, error)
         if(ifg.gt.0) write(6,*) attr_name,'= ',c1buf(1:attrlen)
        endif
        CALL h5aclose_f(attr_id, error)
       enddo
c
      endif
c
      CALL h5eset_auto_f(1,ierr)
c
      return
      end subroutine read_dattr
c=======================================================================
c
c=======================================================================
      SUBROUTINE h5_write(nl,ml,file_id,dsetname,i4buf,r4buf,ifs,
     &                    attr_char0,il,r4attr,ibm, error)
      USE HDF5
c
      implicit none
c
      INTEGER			:: nl,ml
      INTEGER(HID_T)	:: file_id      ! File identifier
      CHARACTER(LEN=80)	:: dsetname		! Dataset name
      integer			:: i4buf(nl*ml)
      REAL				:: r4buf(nl*ml)
      INTEGER			:: ifs
      CHARACTER(LEN=999):: attr_char0	! Attribute data
      INTEGER			:: il
      real(4)			:: r4attr(30)
      integer*2			:: ibm
      integer			:: error		! Error flag
c
      INTEGER(HID_T)	:: dset_id      ! Dataset identifier
      INTEGER(HID_T)	:: dspace_id    ! Dataspace identifier
      INTEGER(HSIZE_T)	:: dims(2)		! Dataset dimensions
      INTEGER			:: rank = 2		! Dataset rank
      CHARACTER(LEN=80)	:: aname        ! Attribute name
      INTEGER(HID_T)	:: attr_id      ! Attribute identifier 
      INTEGER(HID_T)	:: aspace_id    ! Attribute Dataspace identifier 
      INTEGER(HID_T)	:: atype_id     ! Attribute Dataspace identifier 
      INTEGER(HSIZE_T)	:: adims(1)		! Attribute dimension
      INTEGER			:: arank = 1    ! Attribure rank
      INTEGER(SIZE_T)	:: attrlen      ! Length of the attribute string
      CHARACTER(LEN=999):: attr_char	! Attribute data
      INTEGER			:: attr_int		! Attribute data
      REAL				:: attr_real	! Attribute data
      INTEGER(HID_T)	:: plist
      INTEGER(HSIZE_T)	:: chunk(2) =(/256,256/)
      INTEGER			:: i
      INTEGER			:: idmin,idmax,iderr
      REAL				:: rdmin,rdmax,rderr
      INTEGER			:: error0,error1
      logical			:: attr_exists
c
      if(nl.eq.0) then
       nl=ml
       rank=1
      endif
c
      dims(1)=nl
      dims(2)=ml
c
      do while(chunk(1).gt.dims(1))
       chunk(1)=chunk(1)/2
      enddo
      do while(chunk(2).gt.dims(2))
       chunk(2)=chunk(2)/2
      enddo
c
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist, error)
      CALL h5pset_deflate_f(plist, 4, error)
      CALL h5pset_chunk_f(plist, rank, chunk, error)
c
      CALL h5eset_auto_f(0,error0)
      CALL h5dopen_f(file_id, dsetname, dset_id, error1)
      CALL h5eset_auto_f(1,error0)
c
      if(error1.lt.0) then
       if(ifs.eq.2) then
        CALL h5dcreate_f(file_id, dsetname, H5T_STD_U16LE, dspace_id,
     &                  dset_id, error, plist)
       elseif(ifs.eq.-2) then
        CALL h5dcreate_f(file_id, dsetname, H5T_STD_I16LE, dspace_id,
     &                   dset_id, error, plist)
       elseif(ifs.eq.4) then
        CALL h5dcreate_f(file_id, dsetname, H5T_STD_U32LE, dspace_id,
     &                   dset_id, error, plist)
       elseif(ifs.eq.-4) then
        CALL h5dcreate_f(file_id, dsetname, H5T_STD_I32LE, dspace_id,
     &                   dset_id, error, plist)
       elseif(ifs.eq.1) then
        CALL h5dcreate_f(file_id, dsetname, H5T_STD_U8LE, dspace_id,
     &                   dset_id, error, plist)
       elseif(ifs.eq.-1) then
        CALL h5dcreate_f(file_id, dsetname, H5T_STD_I8LE, dspace_id,
     &                   dset_id, error, plist)
       elseif(ifs.eq.0) then
        CALL h5dcreate_f(file_id, dsetname, H5T_IEEE_F32LE, dspace_id,
     &                   dset_id, error, plist)
       endif
      endif
c
      if(ifs.eq.2) then
       idmin=0
       idmax=65534
       iderr=65535
      elseif(ifs.eq.-2) then
       idmin=-32767
       idmax= 32767
       iderr=-32768
      elseif(ifs.eq.4) then
       idmin=0
       idmax=2**16-2+2**16
       iderr=idmax+1
      elseif(ifs.eq.-4) then
       idmax=2**16-1+2**15
       idmin=-idmax
       iderr=idmin-1
      elseif(ifs.eq.1) then
       idmin=0
       idmax=254
       iderr=255
      elseif(ifs.eq.-1) then
       idmin=-127
       idmax= 127
       iderr=-128
      elseif(ifs.eq.0) then
       if(dsetname(16:18).eq.'Lat') then
        rdmin=-90.
        rdmax= 90.
       elseif(dsetname(16:18).eq.'Lon') then
        rdmin=-180.
        rdmax= 180.
       else
        rdmin=-998.
        rdmax=9999.
       endif
       rderr=-999.
      endif
c
      if(ifs.eq.0) then
       CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, r4buf, dims, error)
      else
       CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, i4buf, dims, error)
      endif
c
      do i=1,il
       adims(1)=1
       if     (i==1) then
        aname = "Data_description"
        attr_char=attr_char0
        attrlen=len_trim(attr_char)
        if(dsetname(13:19).eq.'QA_flag') then
         adims(1)=16
         attrlen=attrlen/16
        endif
       elseif (i==2) then
        aname = "Unit"
        if(dsetname(16:23).eq.'Obs_time') then
         attr_char="hour"
        elseif((dsetname(16:21).eq.'Sensor').or.
     &         (dsetname(16:20).eq.'Solar')) then
         attr_char="degree"
        elseif(dsetname(13:21).eq.'Line_msec') then
         attr_char="millisecond"
        elseif(dsetname(13:27).eq.'Land_water_flag') then
         attr_char="NA"
        elseif(dsetname(13:15).eq.'PAR') then
         attr_char="Ein/m^2/day"
        elseif(dsetname(13:15).eq.'BRF') then
         attr_char="NA"
        elseif(dsetname(13:13).eq.'L') then
         attr_char="W/m^2/um/sr"
        elseif(dsetname(13:14).eq.'Tb') then
         attr_char="Kelvin"
        elseif(dsetname(13:16).eq.'CHLA') then
         attr_char="mg/m^3"
        elseif(dsetname(13:16).eq.'NWLR') then
         attr_char="W/m^2/sr/um"
        elseif(dsetname(13:16).eq.'CDOM') then
         attr_char="m^-1"
        elseif(dsetname(13:16).eq.'TSM') then
         attr_char="g m^-3"
        elseif((dsetname(13:15).eq.'aph').or.
     &         (dsetname(13:15).eq.'adg').or.
     &         (dsetname(13:15).eq.'bbp')) then
         attr_char="m^-1"
        else
         attr_char="NA"
        endif
        attrlen=len_trim(attr_char)
       elseif (i==3) then
        aname = "Slope"
        attr_real=r4attr(i)
        attrlen = 0
       elseif (i==4) then
        aname = "Offset"
        attr_real=r4attr(i)
        attrlen = 0
c
       elseif (i==5) then
        if(ifs.eq.0) then
         aname = "Minimum_valid_value"
         attr_real=rdmin
         attrlen =0
        else
         aname = "Minimum_valid_DN"
         attr_int=idmin
         attrlen =-1
        endif
       elseif (i==6) then
        if(ifs.eq.0) then
         aname = "Maximum_valid_value"
         attr_real=rdmax
         attrlen =0
        else
         aname = "Maximum_valid_DN"
         attr_int=idmax
         attrlen =-1
         if(dsetname.eq.'Land_water_flag') attr_int=100
        endif
       elseif (i==7) then
        if(dsetname.eq.'Land_water_flag') then
         aname = "Valid_range"
         attr_char="0(water)-100(land)"
         attrlen=len_trim(attr_char)
        elseif(ifs.eq.0) then
         aname = "Error_value"
         attr_real=rderr
         attrlen =0
        else
         aname = "Error_DN"
         attr_int=iderr
         attrlen =-1
        endif
c
       elseif (i==8) then
        if(dsetname(13:15).eq.'BRF') then
         aname = "Spatial_resolution"
         attr_real=r4attr(i)
         attrlen = 0
        elseif((ifs.eq.-2).or.(ifs.eq.0)) then
         aname = "Resampling_interval"
         attr_int=nint(r4attr(i))
         attrlen =-1
        else
         aname = "Spatial_resolution"
         attr_real=r4attr(i)
         attrlen = 0
        endif
       elseif (i==9) then
        if((ifs.eq.-2).or.(ifs.eq.0)) then
         aname = "Resampling_interval_unit"
         attr_char="pixel"
        elseif(r4attr(8).ge.100.) then
         aname = "Spatial_resolution_unit"
         attr_char="meter"
        elseif(r4attr(8).ge.1) then
         aname = "Resampling_interval_unit"
         attr_char="pixel"
        else
         aname = "Spatial_resolution_unit"
         attr_char="degree"
        endif
        attrlen = len_trim(attr_char)
c
c      ------------------------------------ followings for observation data
       elseif (i==10) then
        aname = "Mask_for_statistics"
        attr_int=ibm
        attrlen =-2
       elseif (i==11) then
        aname = "Center_wavelength"
        attr_real=r4attr(i)
        attrlen =0
       elseif (i==12) then
        aname = "Center_wavelength_unit"
        attr_char="nm"
        attrlen = len_trim(attr_char)
       elseif (i==13) then
        aname = "Band_width"
        attr_real=r4attr(i)
        attrlen =0
       elseif (i==14) then
        aname = "Band_width_unit"
        attr_char="nm"
        attrlen = len_trim(attr_char)
       elseif (i==15) then
        aname = "Band_weighted_TOA_solar_irradiance"
        attr_real=r4attr(i)
        attrlen =0
       elseif (i==16) then
        aname = "Band_weighted_TOA_solar_irradiance_unit"
        attr_char="W/m^2/um"
        attrlen = len_trim(attr_char)
c
       elseif (i==17) then
        if(dsetname(13:16).eq.'NWLR') then
         aname = "Rrs_slope"
        else
         aname = "Slope_reflectance"
        endif
        attr_real=r4attr(i)
        attrlen = 0 
       elseif (i==18) then
        if(dsetname(13:16).eq.'NWLR') then
         aname = "Rrs_offset"
        else
         aname = "Offset_reflectance"
        endif
        attr_real=r4attr(i)
        attrlen = 0 
       elseif (i==19) then
        if(dsetname(13:16).eq.'NWLR') then
         aname = "Rrs_unit"
         attr_char="sr^-1"
        else
         aname = "reflectance_unit"
         attr_char="NA"
        endif
        attrlen = len_trim(attr_char)
       endif
c
       CALL h5eset_auto_f(0,error0)
       CALL h5aexists_f(dset_id, aname, attr_exists, error)
       if( attr_exists ) then
        CALL h5adelete_f(dset_id, aname, error)
       endif
       CALL h5eset_auto_f(1,error0)
c
       CALL h5screate_simple_f(arank, adims, aspace_id, error)
       if(attrlen.gt.0) then
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        CALL h5tset_size_f(atype_id, attrlen, error)
        CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_char, adims, error)
       elseif(attrlen.eq.-2) then
        CALL h5tcopy_f(H5T_STD_U16LE, atype_id, error)
        CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_int, adims, error)
       elseif(attrlen.eq.-1) then
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
        CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_int, adims, error)
       elseif(attrlen.eq.0) then
        CALL h5tcopy_f(H5T_NATIVE_REAL, atype_id, error)
        CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_real, adims, error)
       endif
       CALL h5aclose_f(attr_id, error)
       CALL h5sclose_f(aspace_id, error)
      enddo
c
      CALL h5pclose_f(plist, error)
      CALL h5dclose_f(dset_id, error)
      CALL h5sclose_f(dspace_id, error)
c
      return
      end
c=======================================================================
c
c=======================================================================
      SUBROUTINE h5_write_r8(file_id,dsetname,ml,r8buf, error)
      USE HDF5
c
      implicit none
c
      INTEGER(HID_T)	:: file_id      ! File identifier
      CHARACTER(LEN=80)	:: dsetname		! Dataset name
      INTEGER			:: ml
      REAL(8)			:: r8buf(ml)
      integer			:: error		! Error flag
c
      INTEGER(HID_T)	:: dset_id      ! Dataset identifier
      INTEGER(HID_T)	:: dspace_id    ! Dataspace identifier
      INTEGER			:: rank = 1		! Dataset rank
      INTEGER(HSIZE_T)	:: dims(1)		! Dataset dimensions
      INTEGER(HSIZE_T)	:: chunk(1)
      CHARACTER(LEN=80)	:: aname        ! Attribute name
      INTEGER(HID_T)	:: attr_id      ! Attribute identifier 
      INTEGER(HID_T)	:: aspace_id    ! Attribute Dataspace identifier 
      INTEGER(HID_T)	:: atype_id     ! Attribute Dataspace identifier 
      INTEGER(HSIZE_T)	:: adims(1)		! Attribute dimension
      INTEGER			:: arank = 1    ! Attribure rank
      INTEGER(SIZE_T)	:: attrlen      ! Length of the attribute string
      CHARACTER(LEN=999):: attr_char	! Attribute data
      REAL(4)			:: attr_double	! Attribute data
      INTEGER(HID_T)	:: plist
      INTEGER			:: i,il,error0,error1
      logical			:: attr_exists
c
      dims(1)=ml
      chunk(1)=256
      do i=1,rank
       if(chunk(i).gt.dims(i)) chunk(i)=dims(i)
      enddo
      il=6
c
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist, error)
      CALL h5pset_deflate_f(plist, 4, error)
      CALL h5pset_chunk_f(plist, rank, chunk, error)
c
      CALL h5eset_auto_f(0,error0)
      CALL h5dopen_f(file_id, dsetname, dset_id, error1)
      CALL h5eset_auto_f(1,error0)
c
      if(error1.lt.0) then
       CALL h5dcreate_f(file_id, dsetname, H5T_IEEE_F64LE, dspace_id,
     &                  dset_id, error, plist)
      endif
      CALL h5pclose_f(plist, error)
      CALL h5sclose_f(dspace_id, error)
c
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, r8buf, dims, error)
c
      do i=1,il
       if     (i==1) then
        aname = "Data_description"
        attr_char='TAI93 at each line'
        attrlen=len_trim(attr_char)
       elseif (i==2) then
        aname = "Dim0"
        attr_char="L1B-lines"
        attrlen=len_trim(attr_char)
       elseif (i==3) then
        aname = "Error_value"
        attr_double=-1
        attrlen=-8
       elseif (i==4) then
        aname = "Minimum_valid_DN"
        attr_double=0
        attrlen=-8
       elseif (i==3) then
        aname = "Maximum_valid_value"
        attr_double=999999999
        attrlen=-8
       elseif (i==6) then
        aname = "Unit"
        attr_char="second"
        attrlen=len_trim(attr_char)
       endif
c
       CALL h5eset_auto_f(0,error0)
       CALL h5aexists_f(dset_id, aname, attr_exists, error)
       if( attr_exists ) then
        CALL h5adelete_f(dset_id, aname, error)
       endif
       CALL h5eset_auto_f(1,error0)
c
       adims(1)=1
       CALL h5screate_simple_f(arank, adims, aspace_id, error)
       if(attrlen.gt.0) then
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        CALL h5tset_size_f(atype_id, attrlen, error)
        CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_char, adims, error)
       elseif(attrlen.eq.-8) then
        CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)
        CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_double, adims, error)
       endif
       CALL h5aclose_f(attr_id, error)
       CALL h5sclose_f(aspace_id, error)
      enddo
c
      CALL h5dclose_f(dset_id, error)
c
      return
      end
c=======================================================================
c
c=======================================================================
c     Create Global attributes in "Global_attributes"
      SUBROUTINE gl_attr(file_id,il,filen,stt,ent,ctt,avsn,pvsn,norb,
     &                   error)
c
      USE HDF5
c
      implicit none
c
      INTEGER(HID_T)     :: file_id		! File identifier
      INTEGER			 :: il
      CHARACTER(LEN=50)	 :: filen		! File name
      character(LEN=24)	 :: stt,ent,ctt
      CHARACTER(LEN=4)	 :: avsn		! algorithm version (0~z)
      CHARACTER(LEN=6)	 :: pvsn		! parameter version (0000~zzzz)
      INTEGER			 :: norb		! orbit number
      INTEGER			 :: error		! Error flag
c
      CHARACTER(LEN=4)	 :: fvsn	! product version (0000~zzzz)
      CHARACTER(LEN= 80) :: aname		! Attribute name
      CHARACTER(LEN=999) :: attr_char	! Attribute data
      INTEGER			 :: attr_int	! Attribute data
      CHARACTER(LEN= 80) :: gname		! Group name
      INTEGER(HID_T)	 :: group_id		! File identifier
      INTEGER(SIZE_T)	 :: attrlen		! Length of the attribute string
      INTEGER(HID_T)	 :: attr_id		! Attribute identifier 
      INTEGER			 :: arank = 1	! Attribure rank
      INTEGER(HSIZE_T)	 :: adims(1)= (/1/)	! Attribute dimension
      INTEGER(HID_T)	 :: aspace_id	! Attribute Dataspace identifier 
      INTEGER(HID_T)	 :: atype_id		! Attribute Dataspace identifier 
c
      INTEGER			 :: i
c
      fvsn=avsn(1:1)//pvsn(1:3)
c
      gname='/Global_attributes'
      CALL h5gcreate_f(file_id, gname, group_id, error)
c
      do i=1,il
       if     (i==1) then
        aname = 'Satellite'
        attr_char=
     &  'Global Change Observation Mission - Climate (GCOM-C)'
        attrlen=len_trim(attr_char)
       elseif (i==2) then
        aname = 'Sensor'
        attr_char='Second-generation Global Imager (SGLI)'
        attrlen=len_trim(attr_char)
       elseif (i==3) then
        aname = 'Product_level'
        attr_char='Level-2'
        attrlen=len_trim(attr_char)
       elseif (i==4) then
        aname = 'Product_name'
        attr_char='Atmospheric correction and ocean color'
        attrlen=len_trim(attr_char)
       elseif (i==5) then
        aname = 'Product_version'
        attr_char=fvsn
        attrlen=len_trim(attr_char)
       elseif (i==6) then
        aname = 'Algorithm_version'
        attr_char=avsn
        attrlen=len_trim(attr_char)
       elseif (i==7) then
        aname = 'Parameter_version'
        attr_char=pvsn
        attrlen=len_trim(attr_char)
       elseif (i==8) then
        aname = 'Algorithm_developer'
        attr_char='Japan Aerospace Exploration Agency (JAXA)'
        attrlen=len_trim(attr_char)
       elseif (i==9) then
        aname = 'Dataset_description'
        attr_char='Normalized remote sensing reflectance and Photosynthe
     &tically Available Radiation'
        attrlen=len_trim(attr_char)
       elseif (i==10) then
        aname = 'Product_file_name'
        attr_char=filen
        attrlen=len_trim(attr_char)
       elseif (i==11) then
        if(norb.eq.-999) then
         aname = 'Image_start_time'
        else
         aname = 'Scene_start_time'
        endif
        attr_char=stt(1:21)
        attrlen=len_trim(attr_char)
       elseif (i==12) then
        if(norb.eq.-999) then
         aname = 'Image_end_time'
        else
         aname = 'Scene_end_time'
        endif
        attr_char=ent(1:21)
        attrlen=len_trim(attr_char)
       elseif (i==13) then
        aname = 'Scene_center_time'
        attr_char=ctt(1:21)
        attrlen=len_trim(attr_char)
       elseif (i==14) then
        aname = 'RSP_path_number'
        read(filen(21:23),'(i3)') attr_int
       elseif (i==15) then
        aname = 'Scene_number'
        read(filen(24:25),'(i3)') attr_int
       elseif (i==16) then
        aname = 'Total_orbit_number'
        attr_int=norb
       endif
       CALL h5screate_simple_f(arank, adims, aspace_id, error)
       if(i.le.13) then
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        CALL h5tset_size_f(atype_id, attrlen, error)
        CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_char(1:attrlen), adims,
     &                  error)
       else
        CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
        CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_int, adims, error)
       endif
       CALL h5aclose_f(attr_id, error)
       CALL h5sclose_f(aspace_id, error)
      enddo
c
      CALL h5gclose_f(group_id, error)
c
      return
      end
c=======================================================================
c
c=======================================================================
c     Create Specific attribute groups
      SUBROUTINE sp_attr(file_id,nl,ml,dd,nlg,mlg,ddg,ptime,area,
     &                   ninlst,inlst,qidx,ifanc,ifozn, error)
c
      USE HDF5
c
      implicit none
c
      INTEGER(HID_T)	 :: file_id		! File identifier
      integer 	:: ml,nl
      real(4)  	:: dd
      integer 	:: mlg,nlg
      real(4)  	:: ddg
      character :: ptime*17
      real(4)	:: area(8)	! lat,lon * ul,ur,ll,lr
      integer	:: ninlst
      character :: inlst(501)*48
      character :: inlst44(501)*44
      character :: qidx*4		! Good/Fair/Poor
      integer	:: ifanc,ifozn,error
c
      CHARACTER(LEN= 80) :: gname		! Group name
c      INTEGER(HID_T)	 :: dset_id		! Dataset identifier
c      INTEGER(HID_T)	 :: dspace_id	! Dataspace identifier
      INTEGER(HID_T)	 :: group_id	! Dataspace identifier
      CHARACTER(LEN= 80) :: aname		! Attribute name
      CHARACTER(LEN=999) :: attr_char	! Attribute data
      INTEGER			 :: attr_int	! Attribute data
      REAL				 :: attr_real	! Attribute data
      INTEGER(SIZE_T)	 :: attrlen		! Length of the attribute string
      INTEGER			 :: arank = 1	! Attribure rank
      INTEGER(HSIZE_T)	 :: adims(1)	! Attribute dimension
      INTEGER(HID_T)	 :: aspace_id	! Attribute Dataspace identifier
      INTEGER(HID_T)	 :: atype_id	! Attribute Dataspace identifier  
      INTEGER(HID_T)	 :: attr_id		! Attribute identifier
      integer :: i,k,j,il
c
c     -------------------------------------------- Processing_attributes
      gname="/Processing_attributes"
      CALL h5gcreate_f(file_id, gname, group_id, error)
      il=7
      do i=1,il
       if (i==1) then
        aname = "Processing_organization"
        attr_char="JAXA/GCOM-C science project"
        attrlen=len_trim(attr_char)
        adims(1)=1
       elseif (i==2) then
        aname = "Contact_point"
        attr_char="JAXA/Earth Observation Research Center (EORC)"
        attrlen=len_trim(attr_char)
        adims(1)=1
       elseif (i==3) then
        aname = "Processing_UT"
        attr_char=ptime
        attrlen =17
        adims(1)=1
       elseif (i==4) then
        aname = "Processing_result"
        attr_char=qidx
        attrlen=len_trim(attr_char)
        adims(1)=1
       elseif (i==5) then
        aname = "Input_files"
        attr_char=inlst(1)
        attrlen=len_trim(attr_char)
        adims(1)=ninlst
        if(attrlen.eq.44) then
         do k=1,ninlst
          inlst44(k)=inlst(k)(1:44)
         enddo
        endif
       elseif (i==6) then
        aname = "Ancillary_data_information"
        if    (ifanc.eq.0) then
         attr_char="Standard objective analysis data (JMA GGLA)"
        elseif(ifanc.eq.1) then
         attr_char="Forecast data (JMA GGLF)"
        elseif(ifanc.eq.2) then
         attr_char="Back up ancillary data"
        else
         attr_char="No ancillary data"
        endif
        attrlen=len_trim(attr_char)
        adims(1)=1
       elseif (i==7) then
        aname = "Ozone_data"
        if    (ifozn.eq.0) then
         attr_char="Standard ozone data"
        elseif(ifozn.eq.1) then
         attr_char="Near-real time ozone data"
        elseif(ifozn.eq.2) then
         attr_char="Back-up ozone data"
        elseif(ifozn.eq.3) then
         attr_char="Back-up Near-real time ozone data"
        elseif(ifozn.eq.4) then
         attr_char="2nd back-up ozone data"
        else
         attr_char="No ozone data"
        endif
        attrlen=len_trim(attr_char)
        adims(1)=1
       endif
c
       CALL h5screate_simple_f(arank, adims, aspace_id, error)
       if    (attrlen>0) then		! by character
        CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
        CALL h5tset_size_f(atype_id, attrlen, error)
        CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        if(adims(1).eq.1) then
         CALL h5awrite_f(attr_id, atype_id, attr_char(1:attrlen), adims,
     &                   error)
        elseif(attrlen.eq.44) then
         CALL h5awrite_f(attr_id, atype_id, inlst44, adims, error)
        else
         CALL h5awrite_f(attr_id, atype_id, inlst, adims, error)
        endif
       else					! by real*4
        CALL h5tcopy_f(H5T_NATIVE_REAL, atype_id, error)
        CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                   error)
        CALL h5awrite_f(attr_id, atype_id, attr_real, adims, error)
       endif
       CALL h5aclose_f(attr_id, error)
       CALL h5sclose_f(aspace_id, error)
      enddo
      CALL h5gclose_f(group_id, error)
c
c     -------------------------------------------- Level_1_attributes
      gname="/Level_1_attributes"
      CALL h5gcreate_f(file_id, gname, group_id, error)
      do i=1,3
       adims(1)=1
       if (i==1) then
        aname = "Operation_mode"
        attr_char="OBD"
       elseif (i==2) then
        aname = "Radiometric_calibration"
        attr_char="Original"
       elseif (i==3) then
        aname = "Geometric_calibration"
        attr_char="Original"
       endif
       attrlen=len_trim(attr_char)
       CALL h5screate_simple_f(arank, adims, aspace_id, error)
       CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
       CALL h5tset_size_f(atype_id, attrlen, error)
       CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                  error)
       CALL h5awrite_f(attr_id, atype_id, attr_char(1:attrlen), adims,
     &                 error)
       CALL h5aclose_f(attr_id, error)
       CALL h5sclose_f(aspace_id, error)
      enddo
      CALL h5gclose_f(group_id, error)
c
c     -------------------------------------------- Image & Geometry data
      do j=1,2
       if(j.eq.1) then
        gname="/Image_data"
        il=5
       else
        gname="/Geometry_data"
        il=13
       endif
       CALL h5gcreate_f(file_id, gname, group_id, error)
       do i=1,il
        adims(1)=1
        if (i==1) then
         aname = "Number_of_lines"
         if(j.eq.1) then
          attr_int=ml
         else
          attr_int=mlg
         endif
         attrlen = -1
        elseif (i==2) then
         aname = "Number_of_pixels"
         if(j.eq.1) then
          attr_int=nl
         else
          attr_int=nlg
         endif
         attrlen = -1
        elseif (i==3) then
         aname = "Image_projection"
         if(dd.lt.1.) then
          attr_char=
     &    "EQA (sinusoidal equal area) projection from 0-deg longitude"
          attrlen=len_trim(attr_char)
         else
          attr_char="L1B reference grid"
          attrlen=len_trim(attr_char)
         endif
        elseif (i==4) then
         aname = "Grid_interval"
         if(j.eq.1) then
          attr_real=dd
         else
          attr_real=ddg
         endif
         attrlen = 0
        elseif (i==5) then
         aname = "Grid_interval_unit"
         if(dd.lt.1.) then
          attr_char="deg"
         else
          attr_char="meter"
         endif
         attrlen = len_trim(attr_char)
        elseif (i==6) then
         aname = "Upper_left_longitude"
         attr_real=area(2)
         attrlen = 0
        elseif (i==7) then
         aname = "Upper_left_latitude"
         attr_real=area(1)
         attrlen = 0
        elseif (i==8) then
         aname = "Upper_right_longitude"
         attr_real=area(4)
         attrlen = 0
        elseif (i==9) then
         aname = "Upper_right_latitude"
         attr_real=area(3)
         attrlen = 0
        elseif (i==10) then
         aname = "Lower_left_longitude"
         attr_real=area(6)
         attrlen = 0
        elseif (i==11) then
         aname = "Lower_left_latitude"
         attr_real=area(5)
         attrlen = 0
        elseif (i==12) then
         aname = "Lower_right_longitude"
         attr_real=area(8)
         attrlen = 0
        elseif (i==13) then
         aname = "Lower_right_latitude"
         attr_real=area(7)
         attrlen = 0
        endif
        CALL h5screate_simple_f(arank, adims, aspace_id, error)
        if    (attrlen>0) then		! by character
         CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, error)
         CALL h5tset_size_f(atype_id, attrlen, error)
         CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                    error)
         CALL h5awrite_f(attr_id, atype_id, attr_char, adims, error)
        elseif(attrlen<0) then		! by integer*4
         CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype_id, error)
         CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                    error)
         CALL h5awrite_f(attr_id, atype_id, attr_int, adims, error)
        else							! by real*4
         CALL h5tcopy_f(H5T_NATIVE_REAL, atype_id, error)
         CALL h5acreate_f(group_id, aname, atype_id, aspace_id, attr_id,
     &                    error)
         CALL h5awrite_f(attr_id, atype_id, attr_real, adims, error)
        endif
        CALL h5aclose_f(attr_id, error)
        CALL h5sclose_f(aspace_id, error)
       enddo
       CALL h5gclose_f(group_id, error)
      enddo
c
      return
      end
c=======================================================================
