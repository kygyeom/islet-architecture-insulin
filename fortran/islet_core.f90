!=========================================================================================
!  islet_core.f90  --  CORE islet (beta cells in the centre)
!
!  Islet of Langerhans as a network of coupled Stuart-Landau oscillators (alpha, beta,
!  delta cells) on a 3D HCP lattice, driven by a blood-glucose variable g(t).
!  Source code used for the PLOS ONE manuscript (Core vs Mix islet architecture).
!
!  CORE architecture: beta cells occupy the n_beta lattice sites closest to the
!  islet centre; alpha and delta cells are shuffled at random over the mantle.
!
!  WHAT WAS CHANGED FOR THE PUBLIC RELEASE (nothing else was touched):
!    * English comments were added.
!    * The high-phase input coefficient X, the realization index nt and the random seed
!      are read from the command line instead of being hard-coded / taken from the clock.
!    * The loops over islet size, beta-cell percentage and coupling can optionally be
!      restricted to one value (used for parallel jobs and for the equivalence test).
!    All model equations, constants and the output format are identical to the
!    archived paper-run code.  See README.md, "Equivalence test".
!
!  USAGE
!    islet_core X NT SEED [SIZE] [BRATIO] [EPP]
!      X      high-phase input coefficient: I(t) = X*extamp*g0 for 60 < t <= 120 min.
!             Paper batches: X = 0.2, 0.4, 0.6, 0.8  (nominal G_high = (1+2X)*7 mM =
!             9.8, 12.6, 15.4, 18.2 mM).  Read as default REAL (single precision) on
!             purpose: the archived code used the single-precision literal (e.g. 0.6),
!             so this reproduces the exact same double value after promotion.
!      NT     realization index (1, 2, 3 in the paper). Only used in output file names;
!             different realizations differ through the random numbers.
!      SEED   > 0 : deterministic run. The seed of loop iteration k (k = 1, 2, ...) is
!                   SEED + k - 1. It seeds the intrinsic RANDOM_NUMBER generator
!                   (cell-type placement in init_structure) and sets idum for ran1/gasdev
!                   (natural frequencies and initial oscillator states).
!                   NOTE (archived behaviour, kept): on its very first call ran1 builds
!                   its shuffle table with idum = max(-idum,1) = 1 for any positive seed,
!                   so in the FIRST loop iteration of a program call the oscillator
!                   initial values do not depend on the seed; realizations then differ
!                   through the cell placement. From the second iteration on, ran1
!                   continues from idum = seed of that iteration.
!             <= 0: original behaviour of the paper runs: the seed is taken from the
!                   wall clock (1/100 s resolution).  WARNING: jobs started in the same
!                   1/100 s get identical seeds; use SEED > 0 for parallel runs.
!      SIZE   optional islet-size index sml_size (1..6). 0 or absent = all six.
!      BRATIO optional beta-percentage index b_ratio (1..6). 0 or absent = all six.
!      EPP    optional coupling index epp (1..3 -> K = 0, 0.4, 0.8). 0 or absent = all three.
!
!    Input files read from the working directory: system_parameter.in,
!    structure_parameter.in.  Output: ./DATA_c/  (see "OUTPUT FILES" below).
!
!  LOOP VARIABLES (physical meaning)
!    sml_size = 1..6  islet radius in lattice units  is_size = 1.2*sml_size + 3.8
!                     -> number of cells 743, 1405, 2403, 3791, 5563, 7881
!                     -> diameter label D = 2*5.1um*n^(1/3) = 92, 114, 137, 159, 181, 203 um
!    b_ratio  = 1..6  beta-cell fraction P_beta = 10..60 %  (delta fixed at 10 %,
!                     alpha = 90 % - P_beta)
!    nt       = NT    realization index (independent random initial state / placement)
!    epp      = 1..3  inter-cell (off-diagonal) coupling K = eps = 0.4*(epp-1) = 0, 0.4, 0.8
!    kdd      = 1     self (diagonal) coupling kdia = 0 (fixed in the paper)
!
!  GLUCOSE INPUT  I(t) (variable st), with extamp = 2 and g0 = 7 mM (system_parameter.in)
!       0 <  t <=  60 min : st = 0                    (normal,  G/G0 label 1.0)
!      60 <  t <= 120 min : st = +X  *extamp*g0       (high,    G/G0 label 1 + 2X)
!     120 <  t <= 140 min : st = 0                    (recovery, not analysed)
!     140 <  t <= 200 min : st = -0.2*extamp*g0       (low,     G/G0 label 0.6 = 4.2 mM)
!    The labels are nominal: they are the fixed point of dg/dt = I + g0*ha - g*hb with
!    unit hormone output. The realised g(t) is written to gt_*.dat.
!
!  OUTPUT FILES  (DATA_c/<var>_islet<nt>_<pb>_<sz>_<ep>_<kd>_ext2g0.dat)
!    <pb> = int(10*P_beta) = b_ratio, <sz> = sml_size, <ep> = epp-1 (0,1,2 -> K 0,0.4,0.8),
!    <kd> = kdd-1 = 0.  Two columns (time [min], value), one line every 10 steps
!    (dt = 0.01 min, itmax = 20000 -> 2000 lines, t = 0.1 ... 200.0 min).
!      hb_ : mean beta-cell output per beta cell  hb = (1/N_beta) * sum A_j (1 + cos theta_j)
!            -> the IPA (insulin-production activity) of the paper is the time average
!               of hb over each glucose phase (see scripts/extract_ipa.py).
!      ha_, hd_ : the same for alpha and delta cells;  gt_ : g(t)/g0
!      order_islet*_nonbc_{alpha,beta,delta}_ : Kuramoto order parameter of each type
!    fort.46 (cwd): ha + hb + hd (unit 46 is written without OPEN, legacy behaviour).
!
!  Structure code (init_structure / structure) by Taegeun Song and Danh-Tai Hoang (APCTP);
!  the Monte-Carlo routines inside "structure" are disabled (all calls commented out) and
!  kept only for completeness.
!=========================================================================================
PROGRAM main
implicit none
character(len=32) fname
character(len=32) fname2
character(len=32) fname3, dir, dir2, eeq, kkdi
integer*8 mxdim, nislet, is, info_structure_write
!parameter (mxdim=2500, nislet=10) ! corresponding maximal is_size :   7.5
parameter (mxdim=100000, nislet=1) ! corresponding maximal is_size :   7.5
integer*8 ninfo(nislet,0:4), ncell(nislet), nalpha(nislet), nbeta(nislet), ndelta(nislet)
real*8 is_size(nislet), s_ratio(nislet,3)
real*8 pos(nislet,mxdim,0:4)
integer*8  nj, icase, ii, izref, iz, icount, info, n, i, j, it, itmax, kmax, idum
integer*8  ncell_tmp, nmax
real*8 dt, time, time1, time2, pi, dist
real*8 st, g, dg, eps, extamp, Gam, lamdmax
real*8 rtmp1, rtmp2, rtmp3, rtmp4, rtmp5, rtmp6
real*8 ha, hb, hd, g0, gav, gstdv
real*8 gasdev, ran1, ordertot, rophi1, rophi2, rophi3
real*8 lam(3), pin(3), orderpam(3), ophase(3), kmat(3,3), phi(nislet,3), kdia(3)
real*8 amp(nislet,mxdim), damp(nislet,mxdim), dtheta(nislet,mxdim), theta(nislet,mxdim)
real*8 om(nislet,mxdim), order_nonbc(3), rtmp_nonbc(6)
integer*8 neighbor_info(nislet,mxdim,0:12)
parameter(pi=4.*atan(1.))

integer*8 b_ratio, sml_size, step, nt, ep, kd, ty, epp, kdd
real*8 sml_s(3), nnn(3), test1, eps_p(3), self_p(3), kdi
! --- run-control inputs (public release; read in read_run_args) ---
integer*8 nt_sel, iseed_base, size_lo, size_hi, br_lo, br_hi, ep_lo, ep_hi
real*4    xhigh4    ! X as typed (default REAL, like the archived literal)
real*8    xhigh     ! X promoted to double, used in the input I(t)

!info_structure_write=0       ! write structre
info_structure_write=1       ! do not write structre
!test1 = 0 	! time trace, off 
test1 = 1 	! time trace, on

!call system('rm -r tmp')
!call system('rm -r STRUCTURE')
!call system('mkdir tmp')
!call system('mkdir STRUCTURE')
call read_run_args()          ! X, NT, SEED and optional loop subset
call system('mkdir DATA_c')

!sml_s = [100, 150, 200]
!nnn = [5.5,8.2,10.8]
step = 0
! Archived paper-run loops: sml_size 1..6, b_ratio 1..6, nt 1..3, epp 1..3
! (the archived Core file had 'do epp = 5,5' left over from a test; the paper
!  used epp = 1..3, which is the default here).
! Realizations were run as separate jobs (one NT per program call).
do sml_size = size_lo, size_hi   ! islet size index -> is_size, cell number, D
do b_ratio = br_lo, br_hi        ! beta-cell percentage index -> P_beta = 10*b_ratio %
do nt = nt_sel, nt_sel           ! realization index (file name only)
do epp = ep_lo, ep_hi            ! coupling index -> K = eps = 0.4*(epp-1)
do kdd = 1,1                     ! self coupling index -> kdia = 0.5*(kdd-1) = 0
step = step +1
!-----------------------------------------population ratio & islet size

! Cell-type fractions: alpha = 0.8-0.1*(b_ratio-1), beta = 0.1*b_ratio, delta = 0.1
do is=1, nislet
s_ratio(is,1)=0.8 - 0.1*(b_ratio - 1)
s_ratio(is,2)=0.1 + 0.1*(b_ratio - 1)
s_ratio(is,3)=0.1
enddo

! Islet radius (lattice units); cells are the HCP sites within this radius.
do is=1, nislet
!is_size(is)=6.1
is_size(is)=sml_size*1.2 + 3.8
enddo
!-------------------------------------------------------------------- make structure
call my_timer(time1)
call ran_seed(idum, iseed_base, step)   ! seed = SEED+step-1 (SEED>0) or clock
call get_structure()
call my_timer(time2)
write(6,*) (time2-time1)/60 ,'min elapsed for structure making'
!!------------------------------------------------------------------- parameter setting
! Coupling: off-diagonal (inter-cell) strength K = eps; diagonal kdia = 0
ep = epp - 1
eps = ep*0.4
kd = kdd -1
do i = 1 ,3
kdia(i) = kd*0.5
end do



call my_timer(time1)
call parameter_read()
call neighbor_find()
call my_timer(time2)

kmat=0.
call coupmat(eps,kdia,kmat,info)
pin(1)=Gam ; pin(2)=-Gam ; pin(3)=-Gam



write(6,*) (time2-time1)/60 ,'min elapsed for ready to time trace'
write(6,*) '  '
write(6,*) 'Network ID:', info
write(6,*) 'time window:', itmax*dt
write(6,*) '  '
write(6,*) '   eps:', eps
write(6,*) ' kdia1:', kdia(1)
write(6,*) ' kdia2:', kdia(2)
write(6,*) ' kdia3:', kdia(3)
write(6,*) '    pin:', Gam
write(6,*) 'extamp:', extamp
write(6,*) 'randum_seed', idum

!do is=1, nislet
!do n=1, ncell(is)
!write(6,fmt='(12I5)') (neighbor_info(is,n,i), i=1, neighbor_info(is,n,0))
!enddo
!enddo
!stop
!------------------------------------------------------------------- initial save folder
open(999,file='tmp_fname')
write(999,*) nt 
close(999)
open(999,file='tmp_fname')
read(999,*) fname
close(999)

open(990,file='tmp_nt')
write(990,*) int(10*s_ratio(1,2)) 
close(990)
open(990,file='tmp_nt')
read(990,*) fname2
close(990)

open(990,file='tmp_nt')
write(990,*) int(sml_size)
close(990)
open(990,file='tmp_nt')
read(990,*) fname3
close(990)

open(999,file='tmp_ep')
write(999,*) ep
!write(999,*) int(eps*10)
close(999)
open(999,file='tmp_ep')
read(999,*) eeq
close(999)

open(999,file='tmp_ep')
write(999,*) kd
!write(999,*) int(kdia(1)*10)
close(999)
open(999,file='tmp_ep')
read(999,*) kkdi
close(999) 

! File-name tag: <nt>_<int(10*P_beta)>_<sml_size>  and  _<epp-1>_<kdd-1>
dir = trim(fname)//'_'//trim(fname2)//'_'//trim(fname3)
dir2 = '_'//trim(eeq)//'_'//trim(kkdi)

open(20, file= 'DATA_c'//'/gt_islet'//trim(dir)//trim(dir2)//'_ext2g0.dat')
open(30, file= 'DATA_c'//'/ha_islet'//trim(dir)//trim(dir2)//'_ext2g0.dat')
open(40, file= 'DATA_c'//'/hb_islet'//trim(dir)//trim(dir2)//'_ext2g0.dat')
open(45, file= 'DATA_c'//'/hd_islet'//trim(dir)//trim(dir2)//'_ext2g0.dat')
!open(46, file= 'DATA'//'/total_islet'//trim(dir)//'_ext2g0.dat')
!open(50, file= 'DATA'//'/order_alpha_islet'//trim(dir)//'_ext2g0.dat')
!open(60, file= 'DATA'//'/order_beta_islet'//trim(dir)//'_ext2g0.dat')
!open(70, file= 'DATA'//'/order_delta_islet'//trim(dir)//'_ext2g0.dat')
!open(80, file= 'DATA'//'/order_islet'//trim(dir)//'_ext2g0.dat')
open(90, file= 'DATA_c'//'/order_islet'//trim(dir)//trim(dir2)//'_nonbc_alpha_ext2g0.dat')
open(100,file= 'DATA_c'//'/order_islet'//trim(dir)//trim(dir2)//'_nonbc_beta_ext2g0.dat')
open(110,file= 'DATA_c'//'/order_islet'//trim(dir)//trim(dir2)//'_nonbc_delta_ext2g0.dat')
!------------------------------------------------------------------- initial condition
do is=1, nislet
 do n=1, ncell(is)
   om(is,n)=2.*pi/5.+0.1*gasdev(idum)
 enddo
enddo
if(test1 == 1) then
g=g0+0.1*gasdev(idum) 
do is=1, nislet
 do n=1, ncell(is)
 amp(is,n)=0.5+0.27*ran1(idum)
 theta(is,n)=2.*pi*ran1(idum)
 damp(is,n)=0.2+0.6*ran1(idum)
 dtheta(is,n)=2.*pi*ran1(idum)
 enddo
enddo

!------------------------------------------------------------------------------- time trace
call my_timer(time1)
do it=1, itmax
 time=it*dt
 st = 0
 ! Glucose input I(t): 0 (0-60 min), +X*extamp*g0 (60-120), 0 (120-140),
 ! -0.2*extamp*g0 (140-200).  extamp = 2, g0 = 7 mM.
 if(time.gt.60..and.time.le.120.)  st= xhigh*extamp*g0   ! archived: literal X
 if(time.gt.140..and.time.le.200.) st=-0.2*extamp*g0
 
 do is=1, nislet
  do n=1, ncell(is)
  if(amp(is,n).le.0.)then
  amp(is,n)=-amp(is,n) ; theta(is,n)=theta(is,n)+pi
  endif
  enddo
 enddo

 call lamdacal(lam,g,lamdmax)

 do is=1, nislet
  do n=1, ncell(is)
   if(pos(is,n,0).eq.11.) icase=1                                    !alpha cell 
   if(pos(is,n,0).eq.12.) icase=2                                    ! beta cell 
   if(pos(is,n,0).eq.13.) icase=3                                    !delta cell 

    if(neighbor_info(is,n,0).eq.0)then 
    damp(is,n)=(lam(icase)-abs(amp(is,n))**2)*amp(is,n)
    dtheta(is,n)=om(is,n)+pin(icase)*(g-g0)*cos(theta(is,n))
    endif
 
    if(neighbor_info(is,n,0).ne.0)then 
    rtmp1=0.; rtmp2=0.;
    do j=1, neighbor_info(is,n,0) 
    nj=neighbor_info(is,n,j)
     if(pos(is,nj,0).eq.11.) ii=1                                   
     if(pos(is,nj,0).eq.12.) ii=2                                   
     if(pos(is,nj,0).eq.13.) ii=3                                   
     rtmp1=rtmp1+kmat(ii,icase)*amp(is,nj)*cos(theta(is,nj)-theta(is,n))
     if(amp(is,n).ne.0.) rtmp2=rtmp2+kmat(ii,icase)*amp(is,nj)*sin(theta(is,nj)-theta(is,n))/amp(is,n)
    enddo
    damp(is,n)=(lam(icase)-abs(amp(is,n))**2)*amp(is,n)+rtmp1
    dtheta(is,n)=om(is,n)+pin(icase)*(g-g0)*cos(theta(is,n))+rtmp2
    endif
  enddo
 enddo

!--------------------------------------------------------------Hormon normalized by cell popultion
 ! Mean output per cell of each type: h = (1/N_type) * sum_j A_j (1 + cos theta_j).
 ! hb (beta) is the quantity whose phase averages give the IPA.
 ha=0.;hb=0.; hd=0
 do is=1, nislet
  do n=1, ncell(is)
  if(pos(is,n,0).eq.11.) ha=ha+amp(is,n)*(cos(theta(is,n))+1.)/(nalpha(is)*nislet)
  if(pos(is,n,0).eq.12.) hb=hb+amp(is,n)*(cos(theta(is,n))+1.)/(nbeta(is)*nislet)
  if(pos(is,n,0).eq.13.) hd=hd+amp(is,n)*(cos(theta(is,n))+1.)/(ndelta(is)*nislet)
  enddo
 enddo
!------------------------------------------------------------------------------------------------- 
 ! Glucose balance: dg/dt = g0*ha - g*hb + I(t)
 dg=(ha*g0-hb*g)+st
 call get_order_parameter()
!-------------------------------------------------------write
 ! one line every 10 steps (0.1 min) -> 2000 lines per file
 if(it/10*10.eq.it) write(20,*)  time, g/g0
 !write(20,*)  time, g /g0
 if(it/10*10.eq.it) write(30,*)  time, ha
 if(it/10*10.eq.it) write(40,*)  time, hb
 if(it/10*10.eq.it) write(45,*)  time, hd
 if(it/10*10.eq.it) write(46,*)  time, ha + hb + hd
 !if(it/10*10.eq.it) write(50,*)  time, orderpam(1)
 !if(it/10*10.eq.it) write(60,*)  time, orderpam(2)
 !if(it/10*10.eq.it) write(70,*)  time, orderpam(3)
 !if(it/10*10.eq.it) write(80,*)  time, ordertot
 if(it/10*10.eq.it) write(90,*)   time, order_nonbc(1)
 if(it/10*10.eq.it) write(100,*)  time, order_nonbc(2)
 if(it/10*10.eq.it) write(110,*)  time, order_nonbc(3)
!--------------------------------------------------------update
 g=g+dg*dt
 do is=1, nislet
  do n=1, ncell(is)
  amp(is,n)=amp(is,n)+damp(is,n)*dt
  theta(is,n)=theta(is,n)+dtheta(is,n)*dt
  enddo
 enddo
!--------------------------------------------------------------

enddo
end if

print*, s_ratio(1,2), step ,'loop'


end do ! kdd
end do ! epp

end do ! nt
end do ! b_ratio
end do ! sml_size

call my_timer(time2)
write(6,*) (time2-time1)/60 ,' min elapsed for solving time trace'
print*, '----------------------------'
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '
print*, '      '

contains

subroutine read_run_args()
! Reads:  X NT SEED [SIZE] [BRATIO] [EPP]   (see header).  Public-release addition.
implicit none
character(len=64) arg
integer narg, ios
integer*8 itmp
narg = command_argument_count()
if (narg < 3) then
  write(6,*) 'usage: <exe> X NT SEED [SIZE] [BRATIO] [EPP]'
  stop 1
endif
call get_command_argument(1, arg); read(arg,*,iostat=ios) xhigh4
if (ios /= 0) stop 'bad X'
xhigh = xhigh4                    ! same value as the archived single-precision literal
call get_command_argument(2, arg); read(arg,*,iostat=ios) nt_sel
if (ios /= 0 .or. nt_sel < 1) stop 'bad NT'
call get_command_argument(3, arg); read(arg,*,iostat=ios) iseed_base
if (ios /= 0) stop 'bad SEED'
size_lo = 1; size_hi = 6; br_lo = 1; br_hi = 6; ep_lo = 1; ep_hi = 3
if (narg >= 4) then
  call get_command_argument(4, arg); read(arg,*,iostat=ios) itmp
  if (ios /= 0 .or. itmp < 0 .or. itmp > 6) stop 'bad SIZE'
  if (itmp > 0) then; size_lo = itmp; size_hi = itmp; endif
endif
if (narg >= 5) then
  call get_command_argument(5, arg); read(arg,*,iostat=ios) itmp
  if (ios /= 0 .or. itmp < 0 .or. itmp > 6) stop 'bad BRATIO'
  if (itmp > 0) then; br_lo = itmp; br_hi = itmp; endif
endif
if (narg >= 6) then
  call get_command_argument(6, arg); read(arg,*,iostat=ios) itmp
  if (ios /= 0 .or. itmp < 0 .or. itmp > 3) stop 'bad EPP'
  if (itmp > 0) then; ep_lo = itmp; ep_hi = itmp; endif
endif
write(6,*) 'X =', xhigh, '  NT =', nt_sel, '  SEED =', iseed_base
write(6,*) 'sml_size', size_lo, size_hi, '  b_ratio', br_lo, br_hi, '  epp', ep_lo, ep_hi
return
end subroutine read_run_args

subroutine parameter_read()
implicit none
character(len=64) chtmp
open(900,file='system_parameter.in')
read(900,*) chtmp
read(900,*) chtmp, itmax
read(900,*) chtmp, dt
read(900,*) chtmp, extamp
read(900,*) chtmp
read(900,*) chtmp
read(900,*) chtmp
read(900,*) chtmp!, eps
read(900,*) chtmp!, kdia(1)
read(900,*) chtmp!, kdia(2)
read(900,*) chtmp!, kdia(3)
read(900,*) chtmp
read(900,*) chtmp
read(900,*) chtmp
read(900,*) chtmp, Gam
read(900,*) chtmp
read(900,*) chtmp
read(900,*) chtmp
read(900,*) chtmp, info
read(900,*) chtmp, g0
read(900,*) chtmp, lamdmax
close(900)
return
end subroutine parameter_read

subroutine structure_write()
implicit none
open(999,file='tmp_file')
write(999,fmt='(i8)') is
close(999)
open(999,file='tmp_file')
read(999,*) fname
close(999)
call system('rm tmp_file')
open(100,file='./STRUCTURE/islet'//trim(fname)//'.pdb')
open(200,file='./STRUCTURE/islet'//trim(fname)//'.txt')
do i=1, ncell(is)
if(pos(is,i,0).eq.11)then
write(100,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') 'C',pos(is,i,1),pos(is,i,2),pos(is,i,3),0.
write(200,'(I2,3X,F8.4,3X,F8.4,3X,F8.4)') 12,pos(is,i,1)*4.,pos(is,i,2)*4.,pos(is,i,3)*4.
endif
if(pos(is,i,0).eq.12)then
write(100,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') 'S',pos(is,i,1),pos(is,i,2),pos(is,i,3),0.
write(200,'(I2,3X,F8.4,3X,F8.4,3X,F8.4)') 11,pos(is,i,1)*4.,pos(is,i,2)*4.,pos(is,i,3)*4.
endif
if(pos(is,i,0).eq.13)then
write(100,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') 'Cd',pos(is,i,1),pos(is,i,2),pos(is,i,3),0.
write(200,'(I2,3X,F8.4,3X,F8.4,3X,F8.4)') 13,pos(is,i,1)*4.,pos(is,i,2)*4.,pos(is,i,3)*4.
endif
enddo
close(100)
close(200)
return
end subroutine structure_write
!---------------------------------------------------------find nearlist neighbor
subroutine neighbor_find()
implicit none
neighbor_info=0
do is=1, nislet

 do n=1, ncell(is)
 dist=0.
 icount=0.
 izref=(pos(is,n,3)/sqrt(2./3.)+0.1)/4
 
 do i=1, ncell(is)
 if(i.ne.n)then
  if(pos(is,i,3).eq.pos(is,n,3))then
  dist=sqrt((pos(is,i,1)-pos(is,n,1))**2+(pos(is,i,2)-pos(is,n,2))**2)
   if(dist.lt.5.)then
   icount=icount+1
   neighbor_info(is,n,icount)=i 
   endif
  endif
  if(pos(is,i,3).ne.pos(is,n,3))then
  iz=(pos(is,i,3)/sqrt(2./3.)+0.1)/4
   if(abs(iz-izref).eq.1)then
   dist=sqrt((pos(is,i,1)-pos(is,n,1))**2+(pos(is,i,2)-pos(is,n,2))**2)
   if(dist.lt.4.)then
   icount=icount+1
   neighbor_info(is,n,icount)=i 
   endif
   endif
  endif
 
 endif
 enddo
 neighbor_info(is,n,0)=icount
 enddo
enddo
return
end subroutine neighbor_find
!--------------------------------------------------------- structure make

subroutine get_structure()
implicit none
do is=1, nislet
 call init_structure(is_size, ninfo, pos, mxdim, nislet, is, s_ratio)
 ncell(is)=ninfo(is,0)
 nalpha(is)=ninfo(is,1)
 nbeta(is)=ninfo(is,2)
 ndelta(is)=ninfo(is,3)
 write(6,*) is, 'islet'
 write(6,*) 'total number of cells:', ncell(is) 
 write(6,*) '  number of alpha cells:', nalpha(is), nalpha(is)*1./ncell(is)*1. 
 write(6,*) '  number of beta cells :', nbeta(is), nbeta(is)*1./ncell(is)*1. 
 write(6,*) '  number of delta cells:', ndelta(is), ndelta(is)*1./ncell(is)*1. 
 call structure(mxdim, ninfo, pos, nislet, is)
 if(info_structure_write.eq.0) call structure_write()
enddo ! is
return
end subroutine get_structure
!-----------------------------------------------------------------------------------
!-------------------------------------------------------------------Order parameter
subroutine get_order_parameter()
implicit none
 orderpam=0.; phi=0. ; ordertot=0.

 order_nonbc=0.; rtmp_nonbc=0.

 do is=1, nislet
 rtmp1=0. ; rtmp2=0. ; rtmp3=0. ; rtmp4=0. ; rtmp5=0. ; rtmp6=0.
 rophi1=0.; rophi2=0.; rophi3=0.
  do n=1, ncell(is)
  if(pos(is,n,0).eq.11) rtmp1=rtmp1+cos(theta(is,n))/nalpha(is)
  if(pos(is,n,0).eq.12) rtmp2=rtmp2+cos(theta(is,n))/nbeta(is)
  if(pos(is,n,0).eq.13) rtmp3=rtmp3+cos(theta(is,n))/ndelta(is)
  if(pos(is,n,0).eq.11) rtmp4=rtmp4+sin(theta(is,n))/nalpha(is)
  if(pos(is,n,0).eq.12) rtmp5=rtmp5+sin(theta(is,n))/nbeta(is)
  if(pos(is,n,0).eq.13) rtmp6=rtmp6+sin(theta(is,n))/ndelta(is)

  if(pos(is,n,0).eq.11) rtmp_nonbc(1)=rtmp_nonbc(1)+cos(theta(is,n))/(nislet*nalpha(is))
  if(pos(is,n,0).eq.12) rtmp_nonbc(2)=rtmp_nonbc(2)+cos(theta(is,n))/(nislet*nbeta(is))
  if(pos(is,n,0).eq.13) rtmp_nonbc(3)=rtmp_nonbc(3)+cos(theta(is,n))/(nislet*ndelta(is))
  if(pos(is,n,0).eq.11) rtmp_nonbc(4)=rtmp_nonbc(4)+sin(theta(is,n))/(nislet*nalpha(is))
  if(pos(is,n,0).eq.12) rtmp_nonbc(5)=rtmp_nonbc(5)+sin(theta(is,n))/(nislet*nbeta(is))
  if(pos(is,n,0).eq.13) rtmp_nonbc(6)=rtmp_nonbc(6)+sin(theta(is,n))/(nislet*ndelta(is))

  enddo
 rophi1=sqrt(rtmp1**2+rtmp4**2)
 rophi2=sqrt(rtmp2**2+rtmp5**2)
 rophi3=sqrt(rtmp3**2+rtmp6**2)
 phi(is,1)=phi(is,1)+asin(rtmp4/rophi1)/acos(rtmp1/rophi1)
 phi(is,2)=phi(is,2)+asin(rtmp5/rophi2)/acos(rtmp2/rophi2)
 phi(is,3)=phi(is,3)+asin(rtmp6/rophi3)/acos(rtmp3/rophi3)
 enddo

 order_nonbc(1)=sqrt(rtmp_nonbc(1)**2+rtmp_nonbc(4)**2)
 order_nonbc(2)=sqrt(rtmp_nonbc(2)**2+rtmp_nonbc(5)**2)
 order_nonbc(3)=sqrt(rtmp_nonbc(3)**2+rtmp_nonbc(6)**2)

 rtmp1=0. ; rtmp2=0. ; rtmp3=0. ; rtmp4=0. ; rtmp5=0. ; rtmp6=0.
 do is=1, nislet
  rtmp1=rtmp1+cos(phi(is,1))/nislet
  rtmp2=rtmp2+cos(phi(is,2))/nislet
  rtmp3=rtmp3+cos(phi(is,3))/nislet
  rtmp4=rtmp4+sin(phi(is,1))/nislet
  rtmp5=rtmp5+sin(phi(is,2))/nislet
  rtmp6=rtmp6+sin(phi(is,3))/nislet
 enddo
 orderpam(1)=sqrt(rtmp1**2+rtmp4**2)
 orderpam(2)=sqrt(rtmp2**2+rtmp5**2)
 orderpam(3)=sqrt(rtmp3**2+rtmp6**2)
 ophase(1)=atan(rtmp4/rtmp1)
 ophase(2)=atan(rtmp5/rtmp2)
 ophase(3)=atan(rtmp6/rtmp3)

 rtmp1=0. ; rtmp2=0.
 do i=1, 3
 rtmp1=rtmp1+cos(ophase(i))/3.
 rtmp2=rtmp2+sin(ophase(i))/3.
 enddo
 ordertot=sqrt(rtmp1**2+rtmp2**2)
return
end subroutine get_order_parameter
END PROGRAM main

subroutine coupmat(eps,kdia,kmat,info)
implicit none
integer*8 i, info, itmp1, itmp2, itmp, insign(6)
real*8 eps, kdia(3), kmat(3,3)
insign=0
do i=1, 6
 itmp1=info*(10**(i-6.))
 itmp2=info*(10**(i-7.))
 itmp=itmp1-itmp2*10
 if(mod(itmp,3).eq.1) insign(i)= 1.
 if(mod(itmp,3).eq.2) insign(i)=-1.
 if(mod(itmp,3).eq.0) insign(i)= 0.
enddo
kmat(1,1)=kdia(1)        ; kmat(1,2)=insign(1)*eps  ; kmat(1,3)=insign(3)*eps
kmat(2,1)=insign(2)*eps  ; kmat(2,2)=kdia(2)        ; kmat(2,3)=insign(5)*eps
kmat(3,1)=insign(4)*eps  ; kmat(3,2)=insign(6)*eps  ; kmat(3,3)=kdia(3)
return
end

subroutine lamdacal(lam,g,lamdmax)
implicit none
real*8 lam(3),g,lamdmax
lam(1)= 0.5*(1.+tanh((-g+7.)/5.))   ! lamda alpha
lam(2)= 0.5*(1.+tanh((g-7.)/5.))    ! lamda beta
lam(3)=lamdmax*0.5*(1.+tanh((g-5.)/5.))    ! lamda delta
return
end

subroutine my_timer(ttime)
implicit none
integer*8, parameter :: pd = kind(1.0d0)
real(pd) :: ttime
character(len=8) :: date
character(len=10) :: time
character(len=5) :: zone
integer*8, dimension(8) :: values
integer*8, dimension(12) :: days_in_month
integer*8 :: days
days_in_month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
call date_and_time(date,time,zone,values)
if (mod(values(1),4) == 0) days_in_month(2) = 29
days = sum(days_in_month(1:(values(2)-1)))+values(3)
ttime = 1.0d-3*(values(8)+1.0d3*(values(7)+60*(values(6)+60*(values(5)+24*days))))
return
end

function gasdev(idum)
implicit none
INTEGER*8 idum
REAL*8 gasdev
INTEGER*8 iset
REAL*8 fac,gset,rsq,v1,v2,ran1
SAVE iset,gset
DATA iset/0/
if (iset==0) then
  do
    v1=2.*ran1(idum)-1.
    v2=2.*ran1(idum)-1.
    rsq=v1**2+v2**2
    if(.NOT.(rsq>=1..or.rsq==0.)) EXIT
  end do
  fac=sqrt(-2.*log(rsq)/rsq)
  gset=v1*fac
  gasdev=v2*fac
  iset=1
else
  gasdev=gset
  iset=0
endif
return
end

function ran1(idum)
implicit none
INTEGER*8 idum,IA,IM,IQ,IR,NTAB,NDIV
REAL*8 ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
PARAMETER (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER*8 j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
if (idum.le.0.or.iy.eq.0) then
  idum=max(-idum,1)
!  idum=max(idum,1)
  do 11 j=NTAB+8,1,-1
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    if (j.le.NTAB) iv(j)=idum
11   continue
  iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
END

subroutine ran_seed(idum, iseed_base, step)
! Seeds the intrinsic generator (RANDOM_SEED) and returns idum for ran1/gasdev.
! iseed_base > 0 : deterministic, seed = iseed_base + step - 1 (all other seed
!                  words set to 0 so that the seed is fully defined).
! iseed_base <= 0: archived behaviour, seed from the clock (1/100 s).
implicit none
integer*8 idum, iseed_base, step
integer*8,DIMENSION(8) :: time
integer*4,DIMENSION(50) :: seed
CALL DATE_AND_TIME(values=time)     ! Get the current time
if (iseed_base > 0) then
seed = 0
seed(1) = int(iseed_base + step - 1, 4)
else
seed(1) = time(4)*(360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
endif
CALL RANDOM_SEED(PUT=seed)
idum=seed(1)

return
end



subroutine init_structure(is_size, ninfo, pos, mxdim, nislet, is, s_ratio)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!!! Make HCP structure in 3D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
IMPLICIT NONE
CHARACTER (LEN=3)  :: S
INTEGER*8 ,PARAMETER :: natx=60,naty=60,natz=30
REAL*8    ,PARAMETER :: nul=0.,aaa=4.  
INTEGER*8  :: i,j,k,l,n_cell
REAL*8     :: aaa2,aaa3,x,y,z,x0,y0,z0,r2,r02

integer*8 mxdim, nislet, is, ninfo(nislet,0:4)
real*8 is_size(nislet), pos(nislet,mxdim,0:4)

INTEGER*8  :: icount, n_alpha, n_beta, n_delta
REAL*8     :: ra, rb, rd, s_ratio(nislet,3), temp, ran2
INTEGER*8 , ALLOCATABLE :: ntype(:), idum
! Declare arrays to store distances and ranks
REAL*8, ALLOCATABLE :: distances(:),dx,dy,dz, x_center, y_center, z_center
REAL*8 temp_distance, temp_rank 
INTEGER*8, ALLOCATABLE :: ranks(:), idx


aaa2=aaa*sqrt(3.)/2.
aaa3=aaa*sqrt(8./3.)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      
n_cell=0
r02=is_size(is)

ra=s_ratio(is,1) ; rb=s_ratio(is,2) ; rd=s_ratio(is,3)
!!!!---------------------------------------------------Inserting cell population -> edited Taegeun SONG
icount=0
DO k=1,natz
    DO j=1,naty
        DO i=1,natx
        x0=real(i-1)*aaa
        y0=real(j-1)*aaa2
        z0=real(k-1)*aaa3
        IF (INT(mod(j,2))==1) THEN
            DO l=1,2
            SELECT CASE(l)
            CASE(1)
            x=x0 ; y=y0 ; z=z0
            CASE(2)
            x=x0+aaa/2. ; y=y0+aaa*sqrt(3.)/6. ; z=z0+aaa*sqrt(2./3.) 
            END SELECT
            r2=sqrt((x-aaa*natx/2.)**2.+(y-aaa2*naty/2.)**2.+(z-aaa3*natz/2.)**2.)/aaa
            IF (r2<r02) THEN
                icount=icount+1
            END IF
            ENDDO
        ELSE
            DO l=1,2
            SELECT CASE(l)
            CASE(1)
            x=x0+aaa/2. ; y=y0 ; z=z0                                          
            CASE(2)
            x=x0+aaa   ; y=y0+aaa*sqrt(3.)/6. ; z=z0+aaa*sqrt(2./3.)                                    
            END SELECT
            r2=sqrt((x-aaa*natx/2.)**2.+(y-aaa2*naty/2.)**2.+(z-aaa3*natz/2.)**2.)/aaa
            IF (r2<r02) THEN
                icount=icount+1
            END IF
            ENDDO
        END IF                       
        ENDDO
    ENDDO
  ENDDO
 
 ALLOCATE(ntype(icount))
 n_alpha=ra*icount; n_beta=rb*icount; n_delta=rd*icount

 if(n_alpha+n_beta+n_delta.ne.icount)then
 n_beta=n_beta+(icount-n_alpha-n_beta-n_delta)
 endif

do i=1, int(icount-n_beta+1)
 if(i.le.n_alpha) ntype(i)=11
 if(i.gt.n_alpha) ntype(i)=13
enddo

! Shuffle the cell-type array at random (Fisher-Yates)
do i = int(icount-n_beta), 2, -1
 call random_number(ran2)
 j = (ran2* i) + 1
 temp = ntype(i)
 ntype(i) = ntype(j)
 ntype(j) = temp
end do


!!!!---------------------------------------------------Inserting cell population -> edited Taegeun SONG
DO k=1,natz
    DO j=1,naty
        DO i=1,natx
        x0=real(i-1)*aaa
        y0=real(j-1)*aaa2
        z0=real(k-1)*aaa3
!!!!------------------------------------------------------------------------------------------------      
        IF (INT(mod(j,2))==1) THEN

            DO l=1,2
            SELECT CASE(l)
            CASE(1)
            x=x0 ; y=y0 ; z=z0
            CASE(2)
            x=x0+aaa/2. ; y=y0+aaa*sqrt(3.)/6. ; z=z0+aaa*sqrt(2./3.) 
            END SELECT

            r2=sqrt((x-aaa*natx/2.)**2.+(y-aaa2*naty/2.)**2.+(z-aaa3*natz/2.)**2.)/aaa

            IF (r2<r02) THEN
                n_cell=n_cell+1 
                pos(is,n_cell,1)=x
                pos(is,n_cell,2)=y
                pos(is,n_cell,3)=z
            END IF

            ENDDO
!!!!------------------------------------------------------------------------------------------------  
        ELSE
            DO l=1,2
            SELECT CASE(l)

            CASE(1)
            x=x0+aaa/2. ; y=y0 ; z=z0                                          
            CASE(2)
            x=x0+aaa   ; y=y0+aaa*sqrt(3.)/6. ; z=z0+aaa*sqrt(2./3.)                                    
            END SELECT
            r2=sqrt((x-aaa*natx/2.)**2.+(y-aaa2*naty/2.)**2.+(z-aaa3*natz/2.)**2.)/aaa

            IF (r2<r02) THEN
                n_cell=n_cell+1 
                pos(is,n_cell,1)=x
                pos(is,n_cell,2)=y
                pos(is,n_cell,3)=z
            END IF

            ENDDO

        END IF                       
!!!!------------------------------------------------------------------------------------------------  
        ENDDO
    ENDDO
  ENDDO



! CORE: rank the sites by distance from the centre; the n_beta innermost sites
! become beta cells, the remaining sites take the shuffled alpha/delta types.
ALLOCATE(distances(n_cell))
ALLOCATE(ranks(n_cell))

! Calculate distance from the structure's center
x_center = aaa * natx / 2.0
y_center = aaa2 * naty / 2.0
z_center = aaa3 * natz / 2.0

! Calculate distances and rank the cells based on distance
DO i = 1, n_cell
    x = pos(is, i, 1)
    y = pos(is, i, 2)
    z = pos(is, i, 3)


    dx = x - x_center
    dy = y - y_center
    dz = z - z_center


    distances(i) = SQRT(dx**2 + dy**2 + dz**2)
    ranks(i) = i ! Store the initial index for ranking
ENDDO

! Sort the cells based on distance
DO i = 1, n_cell
    DO j = i + 1, n_cell
        IF (distances(j) < distances(i)) THEN
            temp_distance = distances(i)
            distances(i) = distances(j)
            distances(j) = temp_distance
            
            ! Update the ranking index
            temp_rank = ranks(i)
            ranks(i) = ranks(j)
            ranks(j) = temp_rank
        END IF
    ENDDO
ENDDO

! Assign ntype based on ranks and n_beta, n_alpha, n_delta
DO i = 1, n_cell
    idx = ranks(i) ! Get the original index before sorting
    IF (i.le.n_beta) THEN
    	pos(is, idx, 0) = 12
    ELSE IF (i.gt.n_beta.and.i.le.n_alpha+n_beta) THEN
    	pos(is, idx, 0) = ntype(int(i-n_beta))
    ELSE
    	pos(is, idx, 0) = ntype(int(i-n_beta))
    END IF
ENDDO

! Deallocate arrays
DEALLOCATE(distances)
DEALLOCATE(ranks)



ninfo(is,0)=n_cell
ninfo(is,1)=n_alpha
ninfo(is,2)=n_beta
ninfo(is,3)=n_delta
DEALLOCATE(ntype)
END SUBROUTINE init_structure

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!     Hogil Kim Memorial Building #501 POSTECH,
!     San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!     Personal site: http://hoangdanhtai.com
!-----------------------------------------------------------------------------------------!
!!!   2014.07.29: make structrure HCP with delta cell
!!!               S=11--> alpha, S=12-->beta, S=13 ---> delta
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
subroutine structure(mxdim,ninfo,pos,nislet,is)
IMPLICIT NONE
integer*8 mxdim, nislet, is, ninfo(nislet,0:4)
real*8 pos(nislet,mxdim,0:4), pos_tmp(mxdim,0:4)

CHARACTER (LEN=3)  :: SAt
CHARACTER (LEN=150):: name_data,CONFIG_INI,FIND_NN
CHARACTER (LEN=50) :: name
CHARACTER (LEN=15) :: tmp
CHARACTER (256)    :: Ligne21,Ligne23,Ligne24,Ligne25
REAL*8    ,PARAMETER :: nul=0.,aaa=1.0
INTEGER*8 :: i,j,n_cell,n_max,i_n,ic,S,i_Jbb,n_Jbb,i_Jab,n_Jab,number,i0,i0_n,k,i_n2
INTEGER*8 :: i_delN,i_loop,i_loop1,i_loop2,n_equi1,n_equi2,n_average,i_times,n_times
INTEGER*8 :: nn_total,nn_max,na,nb,n_T,i_T,ia,ib,n_min,nn_min,tmp2,j_ic,name_tmp,j_i
REAL*8    :: r02,r2,n_cell_p,delT,T,Tmin,Tmax,r02max
REAL*8    :: energy,E1_old,E2_old,E1_new,E2_new,Naa,Nbb,Nab,H1a,H1b,H2a,H2b
REAL*8    :: Jaa,Jbb,Jab,del_Jbb,Jbbmax,Jbbmin,del_Jab,Jabmax,Jabmin
REAL*8    :: rdn_ic,rdn_mtp,Naa_av,Nbb_av,Nab_av
REAL*8    :: E_av,Nabc_total,Naa0,Nbb0,Nab0,Nabc0_total
REAL*8    :: pi,alpha,alpha0,OA2,OB2,AB2,tmp1,na_p,real_nb
REAL*8    :: Naa2,Nbb2,Nab2,Naa2_av,Nbb2_av,Nab2_av,FNaa,FNbb,FNab,i_times_compt
INTEGER*8 :: n_Jad,n_Jbd,n_Jdd,id,nd,i_Jad,i_Jbd,i_Jdd,i_case
REAL*8    :: Jad,Jbd,Jdd,Jadmax,Jadmin,Jbdmax,Jbdmin,Jddmax,Jddmin
REAL*8    :: H1d,H2d,Nad0,Nbd0,Ndd0,Nad,Nbd,Ndd,Nad2,Nbd2,Ndd2,FNad,FNbd,FNdd
REAL*8    :: Nad_av,Nbd_av,Ndd_av,Nad2_av,Nbd2_av,Ndd2_av
REAL*8    :: del_Jad,del_Jbd,del_Jdd,nd_p,rdn_case,cost,na_nad
INTEGER*8 ,DIMENSION(:),ALLOCATABLE :: nn,spin1,spin2,nn1,name_a,name_b,name_d
INTEGER*8 ,DIMENSION(:,:),ALLOCATABLE :: name_in,a
REAL*8 ,DIMENSION(:),ALLOCATABLE :: Sa,Sb,Sd,x,y,z,Sa0,Sb0,Sd0
REAL*8 ,DIMENSION(:,:),ALLOCATABLE :: cost_tab,r                    
!!!=======================================================================================
!!!=======================================================================================
!CALL system('rm -r config_ini_3D')
!CALL system('mkdir config_ini_3D')
!CALL system('rm -r config_3D')
!CALL system('mkdir config_3D')
!CALL system('rm *.dat*')

CALL ini_rdm_number()
CALL read_input_parameter()

n_cell=ninfo(is,0)

ALLOCATE(x(n_cell),y(n_cell),z(n_cell))
ALLOCATE(Sa(n_cell),Sb(n_cell),Sd(n_cell))
ALLOCATE(nn(n_cell))
ALLOCATE(nn1(0:n_cell))
ALLOCATE(cost_tab(n_Jbb,n_Jab))
ALLOCATE(spin1(n_cell),spin2(n_cell))
ALLOCATE(r(n_cell,0:n_cell))
ALLOCATE(a(n_cell,0:n_cell))
ALLOCATE(Sa0(n_cell),Sb0(n_cell),Sd0(n_cell))

CALL read_data()

ALLOCATE(name_a(na))
ALLOCATE(name_b(nb))
ALLOCATE(name_d(nd))

!----------------------------------------------------
do i=1, n_cell
do j=0, 4
pos_tmp(i,j)=pos(is,i,j)
enddo
enddo
!-----------------------------------------------------
     
IF (n_T==1) THEN
      delT=0.
ELSE
      delT=(Tmax-Tmin)/real(n_T-1)
END IF

IF (n_Jbb==1) THEN
      del_Jbb=0.
ELSE
      del_Jbb=(Jbbmax-Jbbmin)/real(n_Jbb-1)
END IF

IF (n_Jab==1) THEN
      del_Jab=0.
ELSE
      del_Jab=(Jabmax-Jabmin)/real(n_Jab-1)
END IF

IF (n_Jab==1) THEN
      del_Jab=0.
ELSE
      del_Jab=(Jabmax-Jabmin)/real(n_Jab-1)
END IF

IF (n_Jad==1) THEN
      del_Jad=0.
ELSE
      del_Jad=(Jadmax-Jadmin)/real(n_Jad-1)
END IF

IF (n_Jbd==1) THEN
      del_Jbd=0.
ELSE
      del_Jbd=(Jbdmax-Jbdmin)/real(n_Jbd-1)
END IF

IF (n_Jdd==1) THEN
      del_Jdd=0.
ELSE
      del_Jdd=(Jddmax-Jddmin)/real(n_Jdd-1)
END IF

n_cell_p=real(n_cell)+1.
na_p=real(na)+1.
nd_p=real(nd)+1.
real_nb=real(nb)
      
na_nad=real(na)/real(na+nd)
      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
      
!CALL write_config_ini_3D()
!CALL number_neighbor_ini()
!CALL nearest_neighbor()
!CALL number_interaction_ini()
          
      !!! Doi voi i-consider (i0)
IF (FIND_NN=='YES') THEN
!      CALL neighbor_i0()
!      WRITE(6,*)'Finish find NN of i0'
END IF

!OPEN(unit=21,file='average_thermal.dat')
!OPEN(unit=22,file='del_Nabc2_min.dat')
!OPEN(unit=23,file='energy.dat')
!OPEN(unit=24,file='E_times.dat')             
!OPEN(unit=25,file='average_times.dat')

i_delN=0

DO i_T=1,n_T
!WRITE(6,*)'i_T = ', i_T                       
T=Tmin+delT*real(i_T-1)

DO i_Jbb=1,n_Jbb
!      WRITE(6,*)'i_Jbb=',i_Jbb
      Jbb=Jbbmin+del_Jbb*real(i_Jbb-1)
                  
      DO i_Jab=1,n_Jab
      i_delN=i_delN+1
      Jab=Jabmin+del_Jab*real(i_Jab-1)

      DO i_Jad=1,n_Jad
      Jad=Jadmin+del_Jad*real(i_Jad-1)
      DO i_Jbd=1,n_Jbd
      Jbd=Jbdmin+del_Jbd*real(i_Jbd-1)
      DO i_Jdd=1,n_Jdd
      Jdd=Jddmin+del_Jdd*real(i_Jdd-1)
      
!      CALL load_config_ini()
!      CALL read_name_ab() 

!      CALL equi_lattice1()
!      CALL average_thermal()
       
      !CALL value_thermal()
      !WRITE(*,*)'Naa=',Naa,'Nbb=',Nbb,'Nab=',Nab
            
!      CALL write_config_3D()

      cost_tab(i_Jbb,i_Jab)=cost

      END DO
      END DO
      END DO
      END DO
END DO
END DO

!----------------------------------------------------
do i=1, n_cell
do j=0, 4
pos(is,i,j)=pos_tmp(i,j)
enddo
enddo
!-----------------------------------------------------
     
!WRITE(22,'(2I5,5X,F16.3)')MINLOC(cost_tab),MINVAL(cost_tab)
!WRITE(*,*)Jbb,MINLOC(cost_tab),MINVAL(cost_tab)

!CLOSE(21)
!CLOSE(22)
!CLOSE(23)
!CLOSE(24)
!CLOSE(25)
       
CONTAINS


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_lattice1()
      IMPLICIT NONE

      DO i_loop=1,n_equi1
            CALL equi_lattice()
            !CALL value_thermal()
      END DO

      END SUBROUTINE equi_lattice1

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE equi_lattice()
      IMPLICIT NONE
      
      Naa=0. ; Nbb=0. ; Nab=0. ; Nad=0. ; Nbd=0. ; Ndd=0.
      
      DO j_i=1,nb
            i=name_b(j_i)
 
            !!! Choose a-cell or d-cell?
            CALL random_number(rdn_case)              
            IF (rdn_case<na_nad) THEN    
                  i_case=1
            ELSE
                  i_case=2
            END IF                  
            
            SELECT CASE (i_case)            
            !!!------------------------------------------------------------ 
            !!! su dung ic la a-cell bat ky
            CASE(1)
            j_ic=0
            DO WHILE (j_ic==0)           
                  CALL random_number(rdn_ic)
                  j_ic=int(rdn_ic*na_p)
                  ic=name_a(j_ic)
            END DO

            CALL value_H1()
            CALL value_H2()
    
            E1_old=-Jab*H1a-Jbb*H1b-Jbd*H1d
            E2_old=-Jaa*H2a-Jab*H2b-Jad*H2d
                            
                 Sa(i)=1.  ;  Sb(i)=0.
                 Sa(ic)=0. ;  Sb(ic)=1.

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Jaa*H1a-Jab*H1b-Jad*H1d
                 E2_new=-Jab*H2a-Jbb*H2b-Jbd*H2d

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        name_tmp=name_a(j_ic)
                        name_a(j_ic)=name_b(j_i)
                        name_b(j_i)=name_tmp
                 ELSE         
                        Sa(i)=0.  ;  Sb(i)=1.
                        Sa(ic)=1. ;  Sb(ic)=0.

                 END IF
     
            !!!----------------------------------------------------------- 
            !!! su dung ic la d-cell bat ky
            CASE(2)
            j_ic=0 
            DO WHILE (j_ic==0)
                  CALL random_number(rdn_ic)
                  j_ic=int(rdn_ic*nd_p)
                  ic=name_d(j_ic)
            END DO

            CALL value_H1()
            CALL value_H2()
    
            E1_old=-Jab*H1a-Jbb*H1b-Jbd*H1d
            E2_old=-Jad*H2a-Jbd*H2b-Jdd*H2d
                            
                 Sa(i)=0.  ; Sb(i)=0. ;  Sd(i)=1. 
                 Sa(ic)=0. ; Sb(ic)=1. ; Sd(ic)=0.

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Jad*H1a-Jbd*H1b-Jdd*H1d
                 E2_new=-Jab*H2a-Jbb*H2b-Jbd*H2d

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        name_tmp=name_d(j_ic)
                        name_d(j_ic)=name_b(j_i)
                        name_b(j_i)=name_tmp
                 ELSE         
                        Sa(i)=0.  ; Sb(i)=1. ;  Sd(i)=0. 
                        Sa(ic)=0. ; Sb(ic)=0. ; Sd(ic)=1.

                 END IF                               

            END SELECT
            
      END DO
         
      END SUBROUTINE equi_lattice   

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_thermal() : OK?
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_thermal()
      IMPLICIT NONE
          
      energy=0. ;  Naa=0. ; Nbb=0. ; Nab=0. ; Nad=0. ; Nbd=0. ; Ndd=0.

      DO j_i=1,nb
            i=name_b(j_i)
 
            !!! Choose a-cell or d-cell?
            CALL random_number(rdn_case)              
            IF (rdn_case<na_nad) THEN    
                  i_case=1
            ELSE
                  i_case=2
            END IF                  
            
            SELECT CASE (i_case)            
            !!!------------------------------------------------------------ 
            !!! su dung ic la a-cell bat ky
            CASE(1)
            j_ic=0
            DO WHILE (j_ic==0)           
                  CALL random_number(rdn_ic)
                  j_ic=int(rdn_ic*na_p)
                  ic=name_a(j_ic)
            END DO

            CALL value_H1()
            CALL value_H2()
    
            E1_old=-Jab*H1a-Jbb*H1b-Jbd*H1d
            E2_old=-Jaa*H2a-Jab*H2b-Jad*H2d
                            
                 Sa(i)=1.  ;  Sb(i)=0.
                 Sa(ic)=0. ;  Sb(ic)=1.

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Jaa*H1a-Jab*H1b-Jad*H1d
                 E2_new=-Jab*H2a-Jbb*H2b-Jbd*H2d

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        E1_old=E1_new
                        name_tmp=name_a(j_ic)
                        name_a(j_ic)=name_b(j_i)
                        name_b(j_i)=name_tmp
                 ELSE         
                        Sa(i)=0.  ;  Sb(i)=1.
                        Sa(ic)=1. ;  Sb(ic)=0.

                 END IF
     
            !!!----------------------------------------------------------- 
            !!! su dung ic la d-cell bat ky
            CASE(2)
            j_ic=0 
            DO WHILE (j_ic==0)
                  CALL random_number(rdn_ic)
                  j_ic=int(rdn_ic*nd_p)
                  ic=name_d(j_ic)
            END DO

            CALL value_H1()
            CALL value_H2()
    
            E1_old=-Jab*H1a-Jbb*H1b-Jbd*H1d
            E2_old=-Jad*H2a-Jbd*H2b-Jdd*H2d
                            
                 Sa(i)=0.  ; Sb(i)=0. ;  Sd(i)=1. 
                 Sa(ic)=0. ; Sb(ic)=1. ; Sd(ic)=0.

                 CALL value_H1()
                 CALL value_H2()            
                 E1_new=-Jad*H1a-Jbd*H1b-Jdd*H1d
                 E2_new=-Jab*H2a-Jbb*H2b-Jbd*H2d

                 CALL random_number(rdn_mtp)

                 IF (exp((E1_old+E2_old-E1_new-E2_new)/T) > rdn_mtp) THEN
                        E1_old=E1_new
                        name_tmp=name_d(j_ic)
                        name_d(j_ic)=name_b(j_i)
                        name_b(j_i)=name_tmp
                 ELSE         
                        Sa(i)=0.  ; Sb(i)=1. ;  Sd(i)=0. 
                        Sa(ic)=0. ; Sb(ic)=0. ; Sd(ic)=1.

                 END IF                               

            END SELECT
            
           energy=energy+E1_old
            
      END DO
      energy=energy/2./real_nb

      !!! Tinh Naa, Nbb, Nab, Nad, Nbd, Ndd : -----------------------------------
      DO i=1,n_cell
            CALL value_H1()
            Naa=Naa+Sa(i)*H1a ; Nbb=Nbb+Sb(i)*H1b ; Nab=Nab+Sa(i)*H1b+Sb(i)*H1a
            Nad=Nad+Sa(i)*H1d+Sd(i)*H1a ; Nbd=Nbd+Sb(i)*H1d+Sd(i)*H1b ; Ndd=Ndd+Sd(i)*H1d
      END DO   

      Naa=Naa/2. ; Nbb=Nbb/2. ; Nab=Nab/2. ; Nad=Nad/2. ; Nbd=Nbd/2. ; Ndd=Ndd/2.
      Naa2=Naa*Naa ; Nbb2=Nbb*Nbb ; Nab2=Nab*Nab ; Nad2=Nad*Nad ; Nbd2=Nbd*Nbd ; Ndd2=Ndd*Ndd    

      END SUBROUTINE value_thermal


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE average_thermal() : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE average_thermal()
      IMPLICIT NONE

      E_av=0.; Naa_av=0.; Nbb_av=0. ; Nab_av=0. ; Naa2_av=0.; Nbb2_av=0. ; Nab2_av=0.

      Nad_av=0.; Nbd_av=0. ; Ndd_av=0. ; Nad2_av=0.; Nbd2_av=0. ; Ndd2_av=0.
      
      i_times_compt=0.

      DO i_times=1,n_times

            DO i_loop1=1,n_equi2
                  CALL equi_lattice()
            END DO
            
            DO i_loop2=1,n_average
                  CALL value_thermal()        
                 
                  E_av=E_av+energy
                  Naa_av=Naa_av+Naa ; Nbb_av=Nbb_av+Nbb ; Nab_av=Nab_av+Nab
                  Naa2_av=Naa2_av+Naa2 ; Nbb2_av=Nbb2_av+Nbb2 ; Nab2_av=Nab2_av+Nab2

                  Nad_av=Nad_av+Nad ; Nbd_av=Nbd_av+Nbd ; Ndd_av=Ndd_av+Ndd
                  Nad2_av=Nad2_av+Nad2 ; Nbd2_av=Nbd2_av+Nbd2 ; Ndd2_av=Ndd2_av+Ndd2
                 
                  i_times_compt=i_times_compt+1.
                  !WRITE(Ligne25,*)i_times_compt,Jbb,Jab,Naa,Nbb,Nab
                  !WRITE(25,'(a)') trim(Ligne25)

            END DO

!      WRITE(Ligne24,*) i_times,E_av/real(i_times*n_average)
!      WRITE(24,'(a)') trim(Ligne24)
      
!      WRITE(Ligne25,*)i_times,Jbb,Jab,Naa_av/real(i_times*n_average),Nbb_av/real(i_times*n_average),&
!            Nab_av/real(i_times*n_average)
!      WRITE(25,'(a)') trim(Ligne25)

      END DO

      E_av=E_av/real(n_times*n_average)
      
      Naa_av=Naa_av/real(n_times*n_average)
      Nbb_av=Nbb_av/real(n_times*n_average)
      Nab_av=Nab_av/real(n_times*n_average) 
      
      Naa2_av=Naa2_av/real(n_times*n_average)
      Nbb2_av=Nbb2_av/real(n_times*n_average)
      Nab2_av=Nab2_av/real(n_times*n_average)  

      FNaa=Naa2_av-Naa_av**2.
      FNbb=Nbb2_av-Nbb_av**2.
      FNab=Nab2_av-Nab_av**2.
      
      Nad_av=Nad_av/real(n_times*n_average)
      Nbd_av=Nbd_av/real(n_times*n_average)
      Ndd_av=Ndd_av/real(n_times*n_average) 
      
      Nad2_av=Nad2_av/real(n_times*n_average)
      Nbd2_av=Nbd2_av/real(n_times*n_average)
      Ndd2_av=Ndd2_av/real(n_times*n_average)  

      FNad=Nad2_av-Nad_av**2.
      FNbd=Nbd2_av-Nbd_av**2.
      FNdd=Ndd2_av-Ndd_av**2.

      Nabc_total=(Naa_av+Nbb_av+Nab_av+Nad_av+Nbd_av+Ndd_av)
      
      cost=0.5*((Naa0-Naa_av)**2./FNaa+(Nbb0-Nbb_av)**2./FNbb+(Nab0-Nab_av)**2./FNab&
                  +(Nad0-Nad_av)**2./FNad+(Nbd0-Nbd_av)**2./FNbd+(Ndd0-Ndd_av)**2./FNdd)
                 
!      WRITE(Ligne21,*)Jbb,Jab,Jad,Jbd,Jdd,cost,Nabc_total
!      WRITE(21,'(a)') trim(Ligne21)
!
!      WRITE(Ligne23,*) T,E_av,Naa_av,Nbb_av,Nab_av,Nad_av,Nbd_av,Ndd_av,Nabc_total
!      WRITE(23,'(a)') trim(Ligne23)
     
      END SUBROUTINE average_thermal

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE value_H() : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE value_H1()
      IMPLICIT NONE
      
      !!! value H1 
      H1a=0.  ; H1b=0. ; H1d=0.
      DO j=1,nn(i)
            H1a=H1a+Sa(name_in(i,j))
            H1b=H1b+Sb(name_in(i,j))
            H1d=H1d+Sd(name_in(i,j))
      END DO

      END SUBROUTINE value_H1
!!!================================================
      SUBROUTINE value_H2()
      IMPLICIT NONE

      !!! value H2
      H2a=0.  ; H2b=0. ; H2d=0.
      DO j=1,nn(ic)
            H2a=H2a+Sa(name_in(ic,j))
            H2b=H2b+Sb(name_in(ic,j))
            H2d=H2d+Sd(name_in(ic,j))
      END DO

      END SUBROUTINE value_H2

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  Note: called inside 'structure' after the cell types are already placed. No
!  RANDOM_NUMBER call follows it in this program, so it does not affect outputs.
      SUBROUTINE ini_rdm_number()
      IMPLICIT NONE

      INTEGER*8,DIMENSION(8) :: time
      INTEGER*4,DIMENSION(50) :: seed

      CALL DATE_AND_TIME(values=time)     ! Get the current time
      seed(1) = time(4)*(360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
      CALL RANDOM_SEED(PUT=seed)

      END SUBROUTINE ini_rdm_number 
      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_input_parameter()
      IMPLICIT NONE
      CHARACTER (LEN=150) :: tamp
      OPEN(11,file='structure_parameter.in')
 
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(A10))')   tamp, CONFIG_INI
      READ(11, '(A30,(A10))')   tamp, FIND_NN
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(A10))')   tamp, name_data
      READ(11, '(A30,(F12.6))') tamp, r02
      READ(11, '(A30,(F12.6))') tamp, r02max
      READ(11, '(A30,(I5))')    tamp, i0
      READ(11, '(A30,(F12.6))') tamp, alpha0
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jaa
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jbbmin
      READ(11, '(A30,(F12.6))') tamp, Jbbmax
      READ(11, '(A30,(I5))')    tamp, n_Jbb
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jabmin
      READ(11, '(A30,(F12.6))') tamp, Jabmax
      READ(11, '(A30,(I5))')    tamp, n_Jab
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jadmin
      READ(11, '(A30,(F12.6))') tamp, Jadmax
      READ(11, '(A30,(I5))')    tamp, n_Jad
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jbdmin
      READ(11, '(A30,(F12.6))') tamp, Jbdmax
      READ(11, '(A30,(I5))')    tamp, n_Jbd
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Jddmin
      READ(11, '(A30,(F12.6))') tamp, Jddmax
      READ(11, '(A30,(I5))')    tamp, n_Jdd
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F12.6))') tamp, Tmin
      READ(11, '(A30,(F12.6))') tamp, Tmax
      READ(11, '(A30,(I5))')    tamp, n_T
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I12))')   tamp,n_equi1
      READ(11, '(A30,(I12))')   tamp,n_equi2
      READ(11, '(A30,(I12))')   tamp,n_average
      READ(11, '(A30,(I12))')   tamp,n_times
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 

      END SUBROUTINE read_input_parameter
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Change config as random from data : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE load_config_ini()
      IMPLICIT NONE

 !!!==========================================================
!!! Hinh dang ban dau nhu data goc
      IF (CONFIG_INI=='YES') THEN

      Sa(:)=0. ; Sb(:)=0. ; Sd(:)=0.
      
      DO i=1,n_cell
            Sa(i)=Sa0(i)
            Sb(i)=Sb0(i)
            Sd(i)=Sd0(i)
      END DO

      END IF
      
      END SUBROUTINE load_config_ini
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_ini_3D()
      IMPLICIT NONE

!      OPEN(unit=13,file='config_ini_3D/config_ini_3D.pdb')
      
      DO i=1,n_cell
           
            IF (i==i0) THEN 
                  SAt='H'
!                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
            ELSE          
           
            IF (int(Sb(i))==1) THEN 
                  SAt='S'
!                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
            ELSE
                              
            IF (int(Sa(i))==1) THEN 
                  SAt='C'
!                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
            ELSE
                  SAt='Cd'
!                  WRITE(13,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul

            END IF
            END IF   

            END IF
            
      END DO

      CLOSE(13)

      END SUBROUTINE write_config_ini_3D
      
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE write_config_3D()
      IMPLICIT NONE

!      number=10000000+i_delN
!      WRITE(tmp,'(I8)') number
!      name='_BIO_'//TRIM(tmp)
!      OPEN(unit=14,file='config_3D/config_3D'//trim(name)//'.pdb')
!      OPEN(unit=100,file='config_3D/config_3D'//trim(name)//'.txt')

      DO i=1,n_cell
           
            IF (int(Sb(i))==1) THEN 
                  SAt='S'
!                  WRITE(14,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
!                  WRITE(100,'(I2,3X,F8.4,3X,F8.4,3X,F8.4)')12,x(i)*aaa,y(i)*aaa,z(i)*aaa

!--------------------
pos_tmp(i,0)=12
pos_tmp(i,1)=x(i)*aaa
pos_tmp(i,2)=y(i)*aaa
pos_tmp(i,3)=z(i)*aaa
!--------------------

            ELSE
                              
            IF (int(Sa(i))==1) THEN 
                  SAt='C'
!                  WRITE(14,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
!                  WRITE(100,'(I2,3X,F8.4,3X,F8.4,3X,F8.4)')11,x(i)*aaa,y(i)*aaa,z(i)*aaa
!--------------------
pos_tmp(i,0)=11
pos_tmp(i,1)=x(i)*aaa
pos_tmp(i,2)=y(i)*aaa
pos_tmp(i,3)=z(i)*aaa
!--------------------
            ELSE
                  SAt='Cd'
!                  WRITE(14,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
!                  WRITE(100,'(I2,3X,F8.4,3X,F8.4,3X,F8.4)')13,x(i)*aaa,y(i)*aaa,z(i)*aaa

!--------------------
pos_tmp(i,0)=13
pos_tmp(i,1)=x(i)*aaa
pos_tmp(i,2)=y(i)*aaa
pos_tmp(i,3)=z(i)*aaa
!--------------------

            END IF
            END IF

      END DO

!      CLOSE(14)
!      CLOSE(100)

      END SUBROUTINE write_config_3D

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE Open data and Read_data() Read experimental data : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_data()
      IMPLICIT NONE
  
      Sa(:)=0. ; Sb(:)=0. ; Sd(:)=0.
      na=0 ; nb=0 ; nd=0
      Sa0(:)=0. ; Sb0(:)=0. ; Sd0(:)=0.
      DO i=1,n_cell
!            READ(12,*)S,x(i),y(i),z(i)
      S=pos(is,i,0)
      x(i)=pos(is,i,1)
      y(i)=pos(is,i,2)
      z(i)=pos(is,i,3)
            
            IF ((S==12).or.(S==2)) THEN
                  Sa(i)=0. ; Sb(i)=1. ; Sd(i)=0.
                  Sa0(i)=0. ; Sb0(i)=1. ; Sd0(i)=0.
                  
                  ELSE
                  
                  IF ((S==11).or.(S==5)) THEN
                  Sa(i)=1. ; Sb(i)=0. ;  Sd(i)=0.
                  Sa0(i)=1. ; Sb0(i)=0. ; Sd0(i)=0.
      
                  ELSE
                  IF ((S==13).or.(S==6)) THEN
                  Sa(i)=0. ; Sb(i)=0. ; Sd(i)=1.
                  Sa0(i)=0. ; Sb0(i)=0. ; Sd0(i)=1.

                  END IF
                  END IF
                  
            END IF
          na=na+int(Sa(i)) ; nb=nb+int(Sb(i)) ; nd=nd+int(Sd(i)) 
              
      END DO
      
!      WRITE(6,*)'na,nb,nd=',na,nb,nd

      END SUBROUTINE read_data

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE read_name_ab() : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE read_name_ab()
      IMPLICIT NONE
  
      ia=0 ; ib=0 ; id=0
      name_a(:)=0 ; name_b(:)=0 ; name_d(:)=0

      DO i=1,n_cell
            IF (int(Sa0(i))==1) THEN
                  ia=ia+1
                  name_a(ia)=i 
            ELSE
                  IF (int(Sb0(i))==1) THEN
                  ib=ib+1
                  name_b(ib)=i

                  ELSE

                  IF (int(Sd0(i))==1) THEN
                  id=id+1
                  name_d(id)=i 
                  END IF
 
                  END IF
            END IF
      END DO     

      END SUBROUTINE read_name_ab

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Number nearest neighbor from experimental data : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE number_neighbor_ini()
      IMPLICIT NONE
 
!      OPEN(unit=15,file='nn_i_cell.dat') 
      
      !!!--------------------------------------------------------------
      !!! 
      nn_total=0 ;  nn(:)=0
      
      DO i=1,n_cell
            DO j=1,n_cell
           
            r2=(x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.
           
            IF ((0.1<r2).and.(r2<r02)) THEN
                  nn(i)=nn(i)+1
            END IF

            END DO
                                                !!! nn(i): number neigbors of i cell
            nn_total=nn_total+nn(i)             !!! number neigbors of system
            
!            WRITE(15,*)i,nn(i)

      END DO

      !!!--------------------------------------------------------
      !!! Assign name for nearest neighbor of cell i
      
      ALLOCATE(name_in(n_cell,0:n_cell))
      
      name_in(:,:)=0
      
      DO i=1,n_cell
           i_n=0
           
            DO j=1,n_cell
           
                  r2=(x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.
                 
                  IF ((0.1<r2).and.(r2<r02)) THEN
                        i_n=i_n+1
                        name_in(i,i_n)=j         

                  END IF
                    
            END DO

      END DO

!      CLOSE(15)

      
      END SUBROUTINE number_neighbor_ini

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! Number interaction ini : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE number_interaction_ini()
      IMPLICIT NONE

!      OPEN(unit=16,file='Nab0_ini.dat')    

      nn_total=0
 
      !!!-------------------------------------------------
      !!! maximum of number nearest neighbor
      nn_max=MAXVAL(nn) ; nn_min=MINVAL(nn)
      n_max=0 ; n_min=0
      DO i=1,n_cell
      
            IF (nn(i)==nn_max) THEN                  
                  n_max=n_max+1            
            END IF
            
            IF (nn(i)==nn_min) THEN                  
                  n_min=n_min+1            
            END IF
      
      END DO   
      
      !!!-------------------------------------------------------
      !!! Calculate initial number interaction Naa0, Nbb0, Nab0
      Naa0=0. ; Nbb0=0. ; Nab0=0. ; Nad0=0. ; Nbd0=0. ; Ndd0=0.
      DO i=1,n_cell
            CALL value_H1()
            
            Naa0=Naa0+Sa(i)*H1a ; Nbb0=Nbb0+Sb(i)*H1b ; Nab0=Nab0+Sa(i)*H1b+Sb(i)*H1a
            Nad0=Nad0+Sa(i)*H1d+Sd(i)*H1a ; Nbd0=Nbd0+Sb(i)*H1d+Sd(i)*H1b
            Ndd0=Ndd0+Sd(i)*H1d
      
      END DO 
      Naa0=Naa0/2. ; Nbb0=Nbb0/2. ; Nab0=Nab0/2. ; Nad0=Nad0/2. ; Nbd0=Nbd0/2. ; Ndd0=Ndd0/2.
      Nabc0_total=Naa0+Nbb0+Nab0+Nad0+Nbd0+Ndd0          
      
!      WRITE(16,'(A6,I8,6X,A6,I8,6X,A6,I8,6X,A6,I8)')'na:',na,'nb:',nb,'nd:',nd,'ntotal:',na+nb+nd
!      WRITE(*,'(A6,I8,6X,A6,I8,6X,A6,I8,6X,A6,I8)')'na:',na,'nb:',nb,'nd:',nd,'ntotal:',na+nb+nd
!      
!      WRITE(16,'(A6,F7.2,6X,A6,F7.3,6X,A7,I8,6X,A6,I8,6X,A7,I2,6X,A6,I2)')&
!      'r02:',r02,'nn_av:',2.*real(Nabc0_total)/real(n_cell),'nn_max:',nn_max,'n_max:',n_max,&
!      'nn_min:',nn_min,'n_min:',n_min      
!      WRITE(*,'(A6,F7.2,6X,A6,F7.3,6X,A7,I8,6X,A6,I8,6X,A7,I2,6X,A6,I2)')&
!      'r02:',r02,'nn_av:',2.*real(Nabc0_total)/real(n_cell),'nn_max:',nn_max,'n_max:',n_max,&
!      'nn_min:',nn_min,'n_min:',n_min
!
!      WRITE(16,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1)')&
!               'Naa0:',Naa0,'Nbb0:',Nbb0,'Nab0:',Nab0
!      WRITE(16,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1,6X,A13,F7.1)')&
!               'Ndd0:',Ndd0,'Nad0:',Nad0,'Nbd0:',Nbd0,'N_total:',Nabc0_total
!                        
!      WRITE(*,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1)')&
!               'Naa0:',Naa0,'Nbb0:',Nbb0,'Nab0:',Nab0
!      WRITE(*,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1,6X,A13,F7.1)')&
!               'Ndd0:',Ndd0,'Nad0:',Nad0,'Nbd0:',Nbd0,'N_total:',Nabc0_total
!      
!      CLOSE(16)
!
!      WRITE(*,'(A6,F7.1,6X,A6,F7.1,6X,A6,F7.1)')&
!                'Naa-eq=',Nabc0_total*(real(na)/real(na+nb))**2.,&
!                'Nbb-eq=',Nabc0_total*(real(nb)/real(na+nb))**2.,&
!                'Nab-eq=',Nabc0_total*2.*real(na*nb)/real(na+nb)**2.
!               
     ! WRITE(*,*)'NN i-star:', nn(i0)
      
      
      END SUBROUTINE number_interaction_ini


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! 09.08.2013: nearest neighbor : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE nearest_neighbor()
      IMPLICIT NONE

!      OPEN(unit=31,file='ri_in.dat')
!      OPEN(unit=32,file='config_final_nn.pdb')
!      OPEN(unit=33,file='nn_i_final.dat')

      !!!---------------------------------------------------------
      !!! Tim khoang cach r(i,i_n)
      DO i=1,n_cell
            DO i_n=1,nn(i) 
            r(i,i_n)=(x(i)-x(name_in(i,i_n)))**2.+(y(i)-y(name_in(i,i_n)))**2.+(z(i)-z(name_in(i,i_n)))**2.
            END DO
      END DO

      !!!---------------------------------------------------------
      !!! Sap xep theo thu tu tang dan khoang cach r(i,i_n)
      DO i=1,n_cell
            DO i_n=1,nn(i)-1

            k=i_n

            DO j=i_n+1,nn(i)
                  IF (r(i,j)<r(i,k)) THEN
                        k=j
                  END IF
            END DO

            IF (k/=i_n) THEN
                  tmp1=r(i,i_n)
                  r(i,i_n)=r(i,k)
                  r(i,k)=tmp1

                  tmp2=name_in(i,i_n)
                  name_in(i,i_n)=name_in(i,k)
                  name_in(i,k)=tmp2

            END IF

            END DO
      END DO

      !! Test for i0:
!      DO i_n=1,nn(i0)
!            WRITE(31,*)i_n,r(i0,i_n),name_in(i0,i_n)
!            !WRITE(*,*)i_n,r(i0,i_n),name_in(i0,i_n)
!      END DO

      !!!---------------------------------------------------------
      !!! Kiem tra goc de tim cac nearest neighbor
      pi=acos(-1.)

      a(:,:)=1
      DO i=1,n_cell
      DO i_n=1,nn(i)-1

      IF (a(i,i_n)==1) THEN
      DO i_n2=i_n+1,nn(i)

      OA2=(x(i)-x(name_in(i,i_n)))**2.+(y(i)-y(name_in(i,i_n)))**2.+(z(i)-z(name_in(i,i_n)))**2.
      OB2=(x(i)-x(name_in(i,i_n2)))**2.+(y(i)-y(name_in(i,i_n2)))**2.+(z(i)-z(name_in(i,i_n2)))**2.
      AB2=(x(name_in(i,i_n))-x(name_in(i,i_n2)))**2.+(y(name_in(i,i_n))-y(name_in(i,i_n2)))**2.&
          +(z(name_in(i,i_n))-z(name_in(i,i_n2)))**2.

      alpha=acos((OA2+OB2-AB2)/(2.*sqrt(OA2)*sqrt(OB2)))*180./pi

      !IF ((i==112).and.(i_n==1).and.(i_n2==3)) THEN
       !     WRITE(*,*)'alpha=',alpha
     ! END IF
      IF (alpha<alpha0) THEN
            a(i,i_n2)=0
      END IF

      END DO
      END IF

      END DO
      END DO

      !!!---------------------------------------------------------
      !!! Sap xep lai (theo thu tu khoang cach tang dan) sau khi loai bo
      DO i=1,n_cell
            k=0
            DO i_n=1,nn(i)

            IF (a(i,i_n)==1) THEN
                  k=k+1
                  r(i,k)=r(i,i_n)
                  name_in(i,k)=name_in(i,i_n)
            END IF
            END DO
            
            !! Update nn(i):
            nn(i)=k            

      END DO
      
      !DO i_n=1,nn(i0)
            !WRITE(*,*)'Final:',i_n,r(i0,i_n),name_in(i0,i_n)
      !END DO

      !!!Write for i0

      
           SAt='H'
!           WRITE(32,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!           SAt,x(i0),y(i0),z(i0),nul


          DO i_n=1,nn(i0)
           SAt='Au'
!                  WRITE(32,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(name_in(i0,i_n)),y(name_in(i0,i_n)),z(name_in(i0,i_n)),nul

          END DO
                            
      !!! Number neighbors
      nn1(:)=0

      DO j=1,MAXVAL(nn)
            DO i=1,n_cell
                  IF (nn(i)==j) THEN
                        nn1(j)=nn1(j)+1            
                  END IF
            END DO
!      WRITE(33,*)j,nn1(j)
!      WRITE(*,*)j,nn1(j)

      END DO

      CLOSE(31)
      CLOSE(32)
      CLOSE(33)

      END SUBROUTINE nearest_neighbor

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! neighbor of i0
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE neighbor_i0()
      IMPLICIT NONE
 
      spin1(:)=0 ; spin2(:)=0 ; i0_n=0
       
!      OPEN(unit=17,file='config_ini_i0.pdb')
!      OPEN(unit=18,file='ri0.dat')
           
      DO i=1,n_cell            
            r2=(x(i0)-x(i))**2.+(y(i0)-y(i))**2.+(z(i0)-z(i))**2.
           
            IF ((0.1<r2).and.(r2<r02)) THEN
                  spin1(i)=1
                  i0_n=i0_n+1                
            ELSE                   
            IF ((0.1<r2).and.(r2<r02max)) THEN
                  spin2(i)=1                                                                    
            END IF
            END IF

      END DO
      
!      WRITE(*,*)'Number neigbor of i0 cell (not yet consider angles):',i0_n

!      WRITE(18,*)i0,r02,i0_n

!!!--------------------------------------------------
!!!Write config neigbor of i0      
      DO i=1,n_cell           
            IF (i==i0) THEN 
                  SAt='H'
!                  WRITE(17,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
            END IF            
      END DO
      
      DO i=1,n_cell
            IF (spin1(i)==1) THEN
!                  WRITE(*,*)i

                  SAt='Au'
!                  WRITE(17,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul
            END IF
                              
            IF (spin2(i)==1) THEN 
                  SAt='Cu'
!                  WRITE(17,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
!                        SAt,x(i),y(i),z(i),nul       
            END IF                        
            
      END DO
                              
!      CLOSE(17)
!      CLOSE(18)      

      
      END SUBROUTINE neighbor_i0
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END SUBROUTINE structure
      

       
