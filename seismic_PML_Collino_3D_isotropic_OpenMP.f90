!
! SEISMIC_CPML Version 1.1.1, November 2009.
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               and Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the three-dimensional isotropic elastic wave equation
! using a finite-difference method with classical split Perfectly Matched
! Layer (PML) conditions.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available at the end of this program
! and in file "LICENSE".

  program seismic_PML_Collino_3D_iso

! IMPORTANT : all our CPML codes work fine in single precision as well (which is significantly faster).
!             If you want you can thus force automatic conversion to single precision at compile time
!             or change all the declarations and constants in the code from double precision to single.

  implicit none

!
! 3D explicit PML velocity-stress FD code based upon INRIA report for the 2D case:
!
! Francis Collino and Chrysoula Tsogka
! Application of the PML Absorbing Layer Model to the Linear
! Elastodynamic Problem in Anisotropic Heteregeneous Media
! INRIA Research Report RR-3471, August 1998
! http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
!
! and
!
! @ARTICLE{CoTs01,
! author = {F. Collino and C. Tsogka},
! title = {Application of the {PML} absorbing layer model to the linear elastodynamic
!     problem in anisotropic heterogeneous media},
! journal = {Geophysics},
! year = {2001},
! volume = {66},
! number = {1},
! pages = {294-307}}
!
! PML implemented in the three directions (x, y and z).
!
! Dimitri Komatitsch and Roland Martin, University of Pau, France, April 2007.
!
! The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used.
!
! Parallel implementation based on OpenMP.
! Type for instance "setenv OMP_NUM_THREADS 4" before running in OpenMP if you want 4 tasks.
!
! To display the results as color images in the selected 2D cut plane, use:
!
!   " display image*.gif " or " gimp image*.gif "
!
! or
!
!   " montage -geometry +0+3 -rotate 90 -tile 1x21 image*Vx*.gif allfiles_Vx.gif "
!   " montage -geometry +0+3 -rotate 90 -tile 1x21 image*Vy*.gif allfiles_Vy.gif "
!   then " display allfiles_Vx.gif " or " gimp allfiles_Vx.gif "
!   then " display allfiles_Vy.gif " or " gimp allfiles_Vy.gif "
!

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 101
  integer, parameter :: NY = 641
  integer, parameter :: NZ = 640

! size of a grid cell
  double precision, parameter :: h = 10.d0

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 10

! P-velocity, S-velocity and density
  double precision, parameter :: cp = 3300.d0
  double precision, parameter :: cs = cp / 1.732d0
  double precision, parameter :: rho = 2800.d0
  double precision, parameter :: mu = rho*cs*cs
  double precision, parameter :: lambda = rho*(cp*cp - 2.d0*cs*cs)
  double precision, parameter :: lambda_plus_two_mu = rho*cp*cp

! total number of time steps
  integer, parameter :: NSTEP = 2500

! time step in seconds
  double precision, parameter :: DELTAT = 1.6d-3
  double precision, parameter :: ONE_OVER_DELTAT = 1.d0 / DELTAT

! parameters for the source
  double precision, parameter :: f0 = 7.d0
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.d7

! source
  integer, parameter :: ISOURCE = NX - 2*NPOINTS_PML - 1
  integer, parameter :: JSOURCE = 2 * NY / 3 + 1
  integer, parameter :: KSOURCE = NZ / 2
  double precision, parameter :: xsource = (ISOURCE - 1) * h
  double precision, parameter :: ysource = (JSOURCE - 1) * h
  double precision, parameter :: zsource = (KSOURCE - 1) * h
! angle of source force clockwise with respect to vertical (Y) axis
  double precision, parameter :: ANGLE_FORCE = 135.d0

! receivers
  integer, parameter :: NREC = 2
  double precision, parameter :: xdeb = xsource - 100.d0   ! first receiver x in meters
  double precision, parameter :: ydeb = 2300.d0            ! first receiver y in meters
  double precision, parameter :: zdeb = zsource            ! first receiver y in meters
  double precision, parameter :: xfin = xsource            ! last receiver x in meters
  double precision, parameter :: yfin =  300.d0            ! last receiver y in meters
  double precision, parameter :: zfin =  zsource           ! last receiver y in meters

! display information on the screen from time to time
  integer, parameter :: IT_DISPLAY = 100

! value of PI
  double precision, parameter :: PI = 3.141592653589793238462643d0

! conversion from degrees to radians
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! zero
  double precision, parameter :: ZERO = 0.d0

! large value for maximum
  double precision, parameter :: HUGEVAL = 1.d+30

! velocity threshold above which we consider that the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
  double precision, dimension(NX,NY,NZ) :: vx_1,vx_2,vx_3,&
                                           vy_1,vy_2,vy_3,&
                                           vz_1,vz_2,vz_3,&
                                           sigmaxx_1,sigmaxx_2,sigmaxx_3,&
                                           sigmayy_1,sigmayy_2,sigmayy_3,&
                                           sigmazz_1,sigmazz_2,sigmazz_3,&
                                           sigmaxy_1,sigmaxy_2,&
                                           sigmaxz_1,sigmaxz_3,&
                                           sigmayz_2,sigmayz_3

  double precision, dimension(NX) :: dx_over_two,dx_half_over_two
  double precision, dimension(NY) :: dy_over_two,dy_half_over_two
  double precision, dimension(NZ) :: dz_over_two,dz_half_over_two

! for the source
  double precision a,t,force_x,force_y,force_z,source_term

! for receivers
  double precision xspacerec,yspacerec,zspacerec,distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec,iz_rec
  double precision, dimension(NREC) :: xrec,yrec,zrec
  double precision, dimension(NSTEP,NREC) :: sisvx,sisvy

! for evolution of total energy in the medium
  double precision :: epsilon_xx,epsilon_yy,epsilon_zz,epsilon_xy,epsilon_xz,epsilon_yz
  double precision :: sigmaxx_total,sigmayy_total,sigmazz_total
  double precision :: sigmaxy_total,sigmaxz_total,sigmayz_total
  double precision :: total_energy_kinetic,total_energy_potential
  double precision, dimension(NSTEP) :: total_energy

  integer :: i,j,k,it,irec,iplane

  double precision :: xval,delta,xoriginleft,xoriginright,rcoef,d0,Vsolidnorm,Courant_number,value_dx,value_dy,value_dz,d

! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU

! names of the time stamp files
  character(len=150) outputname

! main I/O file
  integer, parameter :: IOUT = 41

!---
!--- program starts here
!---

!--- define profile of absorption in PML region

! thickness of the layer in meters
  delta = NPOINTS_PML * h

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0 = 3.d0 * cp * log(1.d0/Rcoef) / (2.d0 * delta)

  print *,'d0 = ',d0
  print *

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = delta
  xoriginright = (NX-1)*h - delta

  do i=1,NX

  xval = h*dble(i-1)

  if(xval < xoriginleft) then
    dx_over_two(i) = d0 * ((xoriginleft-xval)/delta)**2
    dx_half_over_two(i) = d0 * ((xoriginleft-xval-h/2.d0)/delta)**2
! fix problem with dx_half_over_two() exactly on the edge
  else if(xval >= 0.9999d0*xoriginright) then
    dx_over_two(i) = d0 * ((xval-xoriginright)/delta)**2
    dx_half_over_two(i) = d0 * ((xval+h/2.d0-xoriginright)/delta)**2
  else
    dx_over_two(i) = 0.d0
    dx_half_over_two(i) = 0.d0
  endif

  enddo

! divide the whole profile by two once and for all
  dx_over_two(:) = dx_over_two(:) / 2.d0
  dx_half_over_two(:) = dx_half_over_two(:) / 2.d0

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = delta
  xoriginright = (NY-1)*h - delta

  do j=1,NY

  xval = h*dble(j-1)

  if(xval < xoriginleft) then
    dy_over_two(j) = d0 * ((xoriginleft-xval)/delta)**2
    dy_half_over_two(j) = d0 * ((xoriginleft-xval-h/2.d0)/delta)**2
! fix problem with dy_half_over_two() exactly on the edge
  else if(xval >= 0.9999d0*xoriginright) then
    dy_over_two(j) = d0 * ((xval-xoriginright)/delta)**2
    dy_half_over_two(j) = d0 * ((xval+h/2.d0-xoriginright)/delta)**2
  else
    dy_over_two(j) = 0.d0
    dy_half_over_two(j) = 0.d0
  endif

  enddo

! divide the whole profile by two once and for all
  dy_over_two(:) = dy_over_two(:) / 2.d0
  dy_half_over_two(:) = dy_half_over_two(:) / 2.d0

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = delta
  xoriginright = (NZ-1)*h - delta

  do k=1,NZ

  xval = h*dble(k-1)

  if(xval < xoriginleft) then
    dz_over_two(k) = d0 * ((xoriginleft-xval)/delta)**2
    dz_half_over_two(k) = d0 * ((xoriginleft-xval-h/2.d0)/delta)**2
! fix problem with dy_half_over_two() exactly on the edge
  else if(xval >= 0.9999d0*xoriginright) then
    dz_over_two(k) = d0 * ((xval-xoriginright)/delta)**2
    dz_half_over_two(k) = d0 * ((xval+h/2.d0-xoriginright)/delta)**2
  else
    dz_over_two(k) = 0.d0
    dz_half_over_two(k) = 0.d0
  endif

  enddo

! divide the whole profile by two once and for all
  dz_over_two(:) = dz_over_two(:) / 2.d0
  dz_half_over_two(:) = dz_half_over_two(:) / 2.d0

! print position of the source
  print *
  print *,'Position of the source:'
  print *
  print *,'x = ',xsource
  print *,'y = ',ysource
  print *,'z = ',zsource
  print *

! define location of receivers
  print *
  print *,'There are ',nrec,' receivers'
  print *
  xspacerec = (xfin-xdeb) / dble(NREC-1)
  yspacerec = (yfin-ydeb) / dble(NREC-1)
  zspacerec = (zfin-zdeb) / dble(NREC-1)
  do irec=1,nrec
    xrec(irec) = xdeb + dble(irec-1)*xspacerec
    yrec(irec) = ydeb + dble(irec-1)*yspacerec
    zrec(irec) = zdeb + dble(irec-1)*zspacerec
  enddo

! find closest grid point for each receiver
  do irec=1,nrec
    dist = HUGEVAL
    do k = 1,NZ
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((h*dble(i-1) - xrec(irec))**2 + (h*dble(j-1) - yrec(irec))**2 + (h*dble(k-1) - zrec(irec))**2)
      if(distval < dist) then
        dist = distval
        ix_rec(irec) = i
        iy_rec(irec) = j
        iz_rec(irec) = k
      endif
    enddo
    enddo
    enddo
    print *,'receiver ',irec,' x_target,y_target,z_target = ',xrec(irec),yrec(irec),zrec(irec)
    print *,'closest grid point found at distance ',dist,' in i,j,k = ',ix_rec(irec),iy_rec(irec),iz_rec(irec)
    print *
  enddo

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
  Courant_number = cp * DELTAT / h
  print *,'Courant number is ',Courant_number
  print *
  if(Courant_number > 1.d0/sqrt(3.d0)) stop 'time step is too large, simulation will be unstable'

! suppress old files (can be commented out if "call system" is missing in your compiler)
! call system('rm -f Vx_*.dat Vy_*.dat Vz_*.dat image*.pnm image*.gif timestamp*')

! initialize arrays
  vx_1(:,:,:) = 0.d0
  vy_1(:,:,:) = 0.d0
  vz_1(:,:,:) = 0.d0

  vx_2(:,:,:) = 0.d0
  vy_2(:,:,:) = 0.d0
  vz_2(:,:,:) = 0.d0

  vx_3(:,:,:) = 0.d0
  vy_3(:,:,:) = 0.d0
  vz_3(:,:,:) = 0.d0

  sigmaxx_1(:,:,:) = 0.d0
  sigmayy_1(:,:,:) = 0.d0
  sigmazz_1(:,:,:) = 0.d0
  sigmaxy_1(:,:,:) = 0.d0
  sigmaxz_1(:,:,:) = 0.d0

  sigmaxx_2(:,:,:) = 0.d0
  sigmayy_2(:,:,:) = 0.d0
  sigmazz_2(:,:,:) = 0.d0
  sigmaxy_2(:,:,:) = 0.d0
  sigmayz_2(:,:,:) = 0.d0

  sigmaxx_3(:,:,:) = 0.d0
  sigmayy_3(:,:,:) = 0.d0
  sigmazz_3(:,:,:) = 0.d0
  sigmaxz_3(:,:,:) = 0.d0
  sigmayz_3(:,:,:) = 0.d0

! initialize seismograms
  sisvx(:,:) = 0.d0
  sisvy(:,:) = 0.d0

! initialize total energy
  total_energy(:) = 0.d0

  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

!---
!---  beginning of time loop
!---

  do it = 1,NSTEP

    print *,'it = ',it

!----------------------
! compute stress sigma
!----------------------

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,d,value_dx,value_dy,value_dz) &
!$OMP SHARED(vx_1,vx_2,vx_3,vy_1,vy_2,vy_3,vz_1,vz_2,vz_3,sigmaxx_1,sigmaxx_2,sigmaxx_3, &
!$OMP sigmayy_1,sigmayy_2,sigmayy_3,sigmazz_1,sigmazz_2,sigmazz_3,dx_half_over_two,dy_over_two,dz_over_two)
do k=2,NZ
  do j = 2,NY
    do i = 1,NX-1

      value_dx = (vx_1(i+1,j,k) - vx_1(i,j,k)) / h &
               + (vx_2(i+1,j,k) - vx_2(i,j,k)) / h &
               + (vx_3(i+1,j,k) - vx_3(i,j,k)) / h

      value_dy = (vy_1(i,j,k) - vy_1(i,j-1,k)) / h &
               + (vy_2(i,j,k) - vy_2(i,j-1,k)) / h &
               + (vy_3(i,j,k) - vy_3(i,j-1,k)) / h

      value_dz = (vz_1(i,j,k) - vz_1(i,j,k-1)) / h &
               + (vz_2(i,j,k) - vz_2(i,j,k-1)) / h &
               + (vz_3(i,j,k) - vz_3(i,j,k-1)) / h

      d = dx_half_over_two(i)

      sigmaxx_1(i,j,k) = ( sigmaxx_1(i,j,k)*(ONE_OVER_DELTAT - d) + lambda_plus_two_mu * value_dx ) / (ONE_OVER_DELTAT + d)

      sigmayy_1(i,j,k) = ( sigmayy_1(i,j,k)*(ONE_OVER_DELTAT - d) + lambda * value_dx ) / (ONE_OVER_DELTAT + d)

      sigmazz_1(i,j,k) = ( sigmazz_1(i,j,k)*(ONE_OVER_DELTAT - d) + lambda * value_dx ) / (ONE_OVER_DELTAT + d)

      d = dy_over_two(j)

      sigmaxx_2(i,j,k) = ( sigmaxx_2(i,j,k)*(ONE_OVER_DELTAT - d) + lambda * value_dy ) / (ONE_OVER_DELTAT + d)

      sigmayy_2(i,j,k) = ( sigmayy_2(i,j,k)*(ONE_OVER_DELTAT - d) + lambda_plus_two_mu * value_dy ) / (ONE_OVER_DELTAT + d)

      sigmazz_2(i,j,k) = ( sigmazz_2(i,j,k)*(ONE_OVER_DELTAT - d) + lambda * value_dy ) / (ONE_OVER_DELTAT + d)

      d = dz_over_two(k)

      sigmaxx_3(i,j,k) = ( sigmaxx_3(i,j,k)*(ONE_OVER_DELTAT - d) + lambda * value_dz ) / (ONE_OVER_DELTAT + d)

      sigmayy_3(i,j,k) = ( sigmayy_3(i,j,k)*(ONE_OVER_DELTAT - d) + lambda * value_dz ) / (ONE_OVER_DELTAT + d)

      sigmazz_3(i,j,k) = ( sigmazz_3(i,j,k)*(ONE_OVER_DELTAT - d) + lambda_plus_two_mu * value_dz ) / (ONE_OVER_DELTAT + d)

    enddo
  enddo
 enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,d,value_dx,value_dy) SHARED(vx_1,vx_2,vx_3,vy_1, &
!$OMP vy_2,vy_3,sigmaxy_1,sigmaxy_2,dy_half_over_two,dx_over_two)
do k=1,NZ
  do j = 1,NY-1
    do i = 2,NX

      value_dx = (vy_1(i,j,k) - vy_1(i-1,j,k)) / h &
               + (vy_2(i,j,k) - vy_2(i-1,j,k)) / h &
               + (vy_3(i,j,k) - vy_3(i-1,j,k)) / h

      value_dy = (vx_1(i,j+1,k) - vx_1(i,j,k)) / h &
               + (vx_2(i,j+1,k) - vx_2(i,j,k)) / h &
               + (vx_3(i,j+1,k) - vx_3(i,j,k)) / h

      d = dx_over_two(i)

      sigmaxy_1(i,j,k) = ( sigmaxy_1(i,j,k)*(ONE_OVER_DELTAT - d) + mu * value_dx ) / (ONE_OVER_DELTAT + d)

      d = dy_half_over_two(j)

      sigmaxy_2(i,j,k) = ( sigmaxy_2(i,j,k)*(ONE_OVER_DELTAT - d) + mu * value_dy ) / (ONE_OVER_DELTAT + d)

    enddo
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,d,value_dx,value_dz) SHARED(vx_1,vx_2,vx_3, &
!$OMP vz_1,vz_2,vz_3,sigmaxz_1,sigmaxz_3,dz_half_over_two,dx_over_two)
do k=1,NZ-1
  do j = 1,NY
    do i = 2,NX

      value_dx = (vz_1(i,j,k) - vz_1(i-1,j,k)) / h &
               + (vz_2(i,j,k) - vz_2(i-1,j,k)) / h &
               + (vz_3(i,j,k) - vz_3(i-1,j,k)) / h

      value_dz = (vx_1(i,j,k+1) - vx_1(i,j,k)) / h &
               + (vx_2(i,j,k+1) - vx_2(i,j,k)) / h &
               + (vx_3(i,j,k+1) - vx_3(i,j,k)) / h

      d = dx_over_two(i)

      sigmaxz_1(i,j,k) = ( sigmaxz_1(i,j,k)*(ONE_OVER_DELTAT - d) + mu * value_dx ) / (ONE_OVER_DELTAT + d)

      d = dz_half_over_two(k)

      sigmaxz_3(i,j,k) = ( sigmaxz_3(i,j,k)*(ONE_OVER_DELTAT - d) + mu * value_dz ) / (ONE_OVER_DELTAT + d)

    enddo
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,d,value_dy,value_dz) SHARED(vy_1,vy_2,vy_3, &
!$OMP vz_1,vz_2,vz_3,sigmayz_2,sigmayz_3,dy_half_over_two,dz_half_over_two)
do k=1,NZ-1
  do j = 1,NY-1
    do i = 1,NX

      value_dy = (vz_1(i,j+1,k) - vz_1(i,j,k)) / h &
               + (vz_2(i,j+1,k) - vz_2(i,j,k)) / h &
               + (vz_3(i,j+1,k) - vz_3(i,j,k)) / h

      value_dz = (vy_1(i,j,k+1) - vy_1(i,j,k)) / h &
               + (vy_2(i,j,k+1) - vy_2(i,j,k)) / h &
               + (vy_3(i,j,k+1) - vy_3(i,j,k)) / h

      d = dy_half_over_two(j)

      sigmayz_2(i,j,k) = ( sigmayz_2(i,j,k)*(ONE_OVER_DELTAT - d) + mu * value_dy ) / (ONE_OVER_DELTAT + d)

      d = dz_half_over_two(k)

      sigmayz_3(i,j,k) = ( sigmayz_3(i,j,k)*(ONE_OVER_DELTAT - d) + mu * value_dz ) / (ONE_OVER_DELTAT + d)

    enddo
  enddo
enddo
!$OMP END PARALLEL DO

!------------------
! compute velocity
!------------------

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,d,value_dx,value_dy,value_dz) SHARED(vx_1,vx_2, &
!$OMP vx_3,sigmaxx_1,sigmaxx_2,sigmaxx_3,sigmaxy_1,sigmaxy_2,sigmaxz_1,sigmaxz_3,dx_over_two,dy_over_two,dz_over_two)
do k = 2,NZ
  do j = 2,NY
    do i = 2,NX

      value_dx = (sigmaxx_1(i,j,k) - sigmaxx_1(i-1,j,k)) / h &
               + (sigmaxx_2(i,j,k) - sigmaxx_2(i-1,j,k)) / h &
               + (sigmaxx_3(i,j,k) - sigmaxx_3(i-1,j,k)) / h

      value_dy = (sigmaxy_1(i,j,k) - sigmaxy_1(i,j-1,k)) / h &
               + (sigmaxy_2(i,j,k) - sigmaxy_2(i,j-1,k)) / h

      value_dz = (sigmaxz_1(i,j,k) - sigmaxz_1(i,j,k-1)) / h &
               + (sigmaxz_3(i,j,k) - sigmaxz_3(i,j,k-1)) / h

      d = dx_over_two(i)

      vx_1(i,j,k) = ( vx_1(i,j,k)*(ONE_OVER_DELTAT - d) + value_dx / rho ) / (ONE_OVER_DELTAT + d)

      d = dy_over_two(j)

      vx_2(i,j,k) = ( vx_2(i,j,k)*(ONE_OVER_DELTAT - d) + value_dy / rho ) / (ONE_OVER_DELTAT + d)

      d = dz_over_two(k)

      vx_3(i,j,k) = ( vx_3(i,j,k)*(ONE_OVER_DELTAT - d) + value_dz / rho ) / (ONE_OVER_DELTAT + d)

    enddo
  enddo
enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,d,value_dx,value_dy,value_dz) SHARED(vy_1,vy_2, &
!$OMP vy_3,sigmayy_1,sigmayy_2,sigmayy_3,sigmaxy_1,sigmaxy_2,sigmayz_2,sigmayz_3,dx_half_over_two,dy_half_over_two,dz_over_two)
do k = 2,NZ
  do j = 1,NY-1
    do i = 1,NX-1

      value_dx = (sigmaxy_1(i+1,j,k) - sigmaxy_1(i,j,k)) / h &
               + (sigmaxy_2(i+1,j,k) - sigmaxy_2(i,j,k)) / h

      value_dy = (sigmayy_1(i,j+1,k) - sigmayy_1(i,j,k)) / h &
               + (sigmayy_2(i,j+1,k) - sigmayy_2(i,j,k)) / h &
               + (sigmayy_3(i,j+1,k) - sigmayy_3(i,j,k)) / h

      value_dz = (sigmayz_2(i,j,k) - sigmayz_2(i,j,k-1)) / h &
               + (sigmayz_3(i,j,k) - sigmayz_3(i,j,k-1)) / h

      d = dx_half_over_two(i)

      vy_1(i,j,k) = ( vy_1(i,j,k)*(ONE_OVER_DELTAT - d) + value_dx / rho ) / (ONE_OVER_DELTAT + d)

      d = dy_half_over_two(j)

      vy_2(i,j,k) = ( vy_2(i,j,k)*(ONE_OVER_DELTAT - d) + value_dy / rho ) / (ONE_OVER_DELTAT + d)

      d = dz_over_two(k)

      vy_3(i,j,k) = ( vy_3(i,j,k)*(ONE_OVER_DELTAT - d) + value_dz / rho ) / (ONE_OVER_DELTAT + d)

    enddo
  enddo
 enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,d,value_dx,value_dy,value_dz) SHARED(vz_1,vz_2, &
!$OMP vz_3,sigmazz_1,sigmazz_2,sigmazz_3,sigmaxz_1,sigmaxz_3,sigmayz_2,sigmayz_3,dx_half_over_two,dy_over_two,dz_half_over_two)
do k = 1,NZ-1
  do j = 2,NY
    do i = 1,NX-1

      value_dx = (sigmaxz_1(i+1,j,k) - sigmaxz_1(i,j,k)) / h &
               + (sigmaxz_3(i+1,j,k) - sigmaxz_3(i,j,k)) / h

      value_dy = (sigmayz_2(i,j,k) - sigmayz_2(i,j-1,k)) / h &
               + (sigmayz_3(i,j,k) - sigmayz_3(i,j-1,k)) / h

      value_dz = (sigmazz_1(i,j,k+1) - sigmazz_1(i,j,k)) / h &
               + (sigmazz_2(i,j,k+1) - sigmazz_2(i,j,k)) / h &
               + (sigmazz_3(i,j,k+1) - sigmazz_3(i,j,k)) / h

      d = dx_half_over_two(i)

      vz_1(i,j,k) = ( vz_1(i,j,k)*(ONE_OVER_DELTAT - d) + value_dx / rho ) / (ONE_OVER_DELTAT + d)

      d = dy_over_two(j)

      vz_2(i,j,k) = ( vz_2(i,j,k)*(ONE_OVER_DELTAT - d) + value_dy / rho ) / (ONE_OVER_DELTAT + d)

      d = dz_half_over_two(k)

      vz_3(i,j,k) = ( vz_3(i,j,k)*(ONE_OVER_DELTAT - d) + value_dz / rho ) / (ONE_OVER_DELTAT + d)

    enddo
  enddo
 enddo
!$OMP END PARALLEL DO

! add the source (force vector located at a given grid point)
  a = pi*pi*f0*f0
  t = dble(it-1)*DELTAT

! Gaussian
! source_term = factor * exp(-a*(t-t0)**2)

! first derivative of a Gaussian
  source_term = - factor * 2.d0*a*(t-t0)*exp(-a*(t-t0)**2)

! Ricker source time function (second derivative of a Gaussian)
! source_term = factor * (1.d0 - 2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)

  force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term
  force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term
  force_z = 0.d0

! define location of the source
  i = ISOURCE
  j = JSOURCE
  k = NZ / 2

! add the source to one of the two components of the split field
  vx_1(i,j,k) = vx_1(i,j,k) + force_x * DELTAT / rho
  vy_1(i,j,k) = vy_1(i,j,k) + force_y * DELTAT / rho

! implement Dirichlet boundary conditions on the six edges of the grid

!$OMP PARALLEL WORKSHARE
! xmin
  vx_1(1,:,:) = 0.d0
  vy_1(1,:,:) = 0.d0
  vz_1(1,:,:) = 0.d0

  vx_2(1,:,:) = 0.d0
  vy_2(1,:,:) = 0.d0
  vz_2(1,:,:) = 0.d0

  vx_3(1,:,:) = 0.d0
  vy_3(1,:,:) = 0.d0
  vz_3(1,:,:) = 0.d0

! xmax
  vx_1(NX,:,:) = 0.d0
  vy_1(NX,:,:) = 0.d0
  vz_1(NX,:,:) = 0.d0

  vx_2(NX,:,:) = 0.d0
  vy_2(NX,:,:) = 0.d0
  vz_2(NX,:,:) = 0.d0

  vx_3(NX,:,:) = 0.d0
  vy_3(NX,:,:) = 0.d0
  vz_3(NX,:,:) = 0.d0

! ymin
  vx_1(:,1,:) = 0.d0
  vy_1(:,1,:) = 0.d0
  vz_1(:,1,:) = 0.d0

  vx_2(:,1,:) = 0.d0
  vy_2(:,1,:) = 0.d0
  vz_2(:,1,:) = 0.d0

  vx_3(:,1,:) = 0.d0
  vy_3(:,1,:) = 0.d0
  vz_3(:,1,:) = 0.d0

! ymax
  vx_1(:,NY,:) = 0.d0
  vy_1(:,NY,:) = 0.d0
  vz_1(:,NY,:) = 0.d0

  vx_2(:,NY,:) = 0.d0
  vy_2(:,NY,:) = 0.d0
  vz_2(:,NY,:) = 0.d0

  vx_3(:,NY,:) = 0.d0
  vy_3(:,NY,:) = 0.d0
  vz_3(:,NY,:) = 0.d0

! zmin
  vx_1(:,:,1) = 0.d0
  vy_1(:,:,1) = 0.d0
  vz_1(:,:,1) = 0.d0

  vx_2(:,:,1) = 0.d0
  vy_2(:,:,1) = 0.d0
  vz_2(:,:,1) = 0.d0

  vx_3(:,:,1) = 0.d0
  vy_3(:,:,1) = 0.d0
  vz_3(:,:,1) = 0.d0

! zmax
  vx_1(:,:,NZ) = 0.d0
  vy_1(:,:,NZ) = 0.d0
  vz_1(:,:,NZ) = 0.d0

  vx_2(:,:,NZ) = 0.d0
  vy_2(:,:,NZ) = 0.d0
  vz_2(:,:,NZ) = 0.d0

  vx_3(:,:,NZ) = 0.d0
  vy_3(:,:,NZ) = 0.d0
  vz_3(:,:,NZ) = 0.d0
!$OMP END PARALLEL WORKSHARE

! store seismograms
  do irec = 1,NREC
    sisvx(it,irec) = vx_1(ix_rec(irec),iy_rec(irec),iz_rec(irec)) + &
      vx_2(ix_rec(irec),iy_rec(irec),iz_rec(irec)) + vx_3(ix_rec(irec),iy_rec(irec),iz_rec(irec))
    sisvy(it,irec) = vy_1(ix_rec(irec),iy_rec(irec),iz_rec(irec)) + &
      vy_2(ix_rec(irec),iy_rec(irec),iz_rec(irec)) + vy_3(ix_rec(irec),iy_rec(irec),iz_rec(irec))
  enddo

! compute total energy in the medium (without the PML layers)

  total_energy_kinetic = ZERO
  total_energy_potential = ZERO

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,sigmaxx_total,sigmayy_total, &
!$OMP sigmazz_total,sigmaxy_total,sigmaxz_total,sigmayz_total,epsilon_xx,epsilon_yy,epsilon_zz,epsilon_xy,epsilon_xz,epsilon_yz) &
!$OMP SHARED(vx_1,vx_2,vx_3,vy_1,vy_2,vy_3,vz_1,vz_2,vz_3,sigmaxx_1,sigmaxx_2, &
!$OMP sigmaxx_3,sigmayy_1,sigmayy_2,sigmayy_3,sigmazz_1,sigmazz_2,sigmazz_3, &
!$OMP sigmaxy_1,sigmaxy_2,sigmaxz_1,sigmaxz_3,sigmayz_2,sigmayz_3) REDUCTION(+:total_energy_kinetic,total_energy_potential)
  do k = NPOINTS_PML+1, NZ-NPOINTS_PML
    do j = NPOINTS_PML+1, NY-NPOINTS_PML
      do i = NPOINTS_PML+1, NX-NPOINTS_PML

! compute kinetic energy first, defined as 1/2 rho ||v||^2
! in principle we should use rho_half_x_half_y instead of rho for vy
! in order to interpolate density at the right location in the staggered grid cell
! but in a homogeneous medium we can safely ignore it
  total_energy_kinetic = total_energy_kinetic + 0.5d0 * rho*( &
      (vx_1(i,j,k) + vx_2(i,j,k) + vx_3(i,j,k))**2 + &
      (vy_1(i,j,k) + vy_2(i,j,k) + vy_3(i,j,k))**2 + &
      (vz_1(i,j,k) + vz_2(i,j,k) + vz_3(i,j,k))**2)

! add potential energy, defined as 1/2 epsilon_ij sigma_ij
! in principle we should interpolate the medium parameters at the right location
! in the staggered grid cell but in a homogeneous medium we can safely ignore it

! compute total field from split components
      sigmaxx_total = sigmaxx_1(i,j,k) + sigmaxx_2(i,j,k) + sigmaxx_3(i,j,k)
      sigmayy_total = sigmayy_1(i,j,k) + sigmayy_2(i,j,k) + sigmayy_3(i,j,k)
      sigmazz_total = sigmazz_1(i,j,k) + sigmazz_2(i,j,k) + sigmazz_3(i,j,k)
      sigmaxy_total = sigmaxy_1(i,j,k) + sigmaxy_2(i,j,k)
      sigmaxz_total = sigmaxz_1(i,j,k) + sigmaxz_3(i,j,k)
      sigmayz_total = sigmayz_2(i,j,k) + sigmayz_3(i,j,k)

      epsilon_xx = (2.d0*(lambda + mu) * sigmaxx_total - lambda * sigmayy_total -lambda*sigmazz_total) / &
               (2.d0 * mu * (3.d0*lambda + 2.d0*mu))
      epsilon_yy = (2.d0*(lambda + mu) * sigmayy_total - lambda * sigmaxx_total -lambda*sigmazz_total) / &
               (2.d0 * mu * (3.d0*lambda + 2.d0*mu))
      epsilon_zz = (2.d0*(lambda + mu) * sigmazz_total - lambda * sigmaxx_total -lambda*sigmayy_total) / &
               (2.d0 * mu * (3.d0*lambda + 2.d0*mu))
      epsilon_xy = sigmaxy_total / (2.d0 * mu)
      epsilon_xz = sigmaxz_total / (2.d0 * mu)
      epsilon_yz = sigmayz_total / (2.d0 * mu)

      total_energy_potential = total_energy_potential + &
        0.5d0 * (epsilon_xx * sigmaxx_total + epsilon_yy * sigmayy_total + &
        epsilon_yy * sigmayy_total+ 2.d0 * epsilon_xy * sigmaxy_total + &
        2.d0*epsilon_xz * sigmaxz_total+2.d0*epsilon_yz * sigmayz_total)

      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

  total_energy(it) = total_energy_kinetic + total_energy_potential

! output information
  if(mod(it,IT_DISPLAY) == 0 .or. it == 5) then

      Vsolidnorm = maxval(sqrt((vx_1 + vx_2 + vx_3)**2 + (vy_1 + vy_2 + vy_3)**2+(vz_1 + vz_2 + vz_3)**2))

      print *,'Time step # ',it
      print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
      print *,'Max norm velocity vector V (m/s) = ',Vsolidnorm
      print *,'Total energy = ',total_energy(it)
! check stability of the code, exit if unstable
      if(Vsolidnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
    iplane=1

! count elapsed wall-clock time
    call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
    time_end = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! elapsed time since beginning of the simulation
    tCPU = time_end - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(*,*) 'Elapsed time in seconds = ',tCPU
    write(*,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(*,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    write(*,*)

! write time stamp file to give information about progression of simulation
    write(outputname,"('timestamp',i6.6)") it
    open(unit=IOUT,file=outputname,status='unknown')
    write(IOUT,*) 'Time step # ',it
    write(IOUT,*) 'Time: ',sngl((it-1)*DELTAT),' seconds'
    write(IOUT,*) 'Max norm velocity vector V (m/s) = ',Vsolidnorm
    write(IOUT,*) 'Total energy = ',total_energy(it)
    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    close(IOUT)

! save seismograms
    print *,'saving seismograms'
    print *
    call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT)

    call create_color_image(vx_1(:,:,NZ/2) + vx_2(:,:,NZ/2) + vx_3(:,:,NZ/2),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,.true.,.true.,.true.,.true.,1)

    call create_color_image(vy_1(:,:,NZ/2) + vy_2(:,:,NZ/2) +vy_3(:,:,NZ/2),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,.true.,.true.,.true.,.true.,2)

    endif

  enddo   ! end of time loop

! save seismograms
  call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT)

! save total energy
  open(unit=20,file='energy.dat',status='unknown')
  do it = 1,NSTEP
    write(20,*) sngl(dble(it-1)*DELTAT),total_energy(it)
  enddo
  close(20)

! create script for Gnuplot for total energy
  open(unit=20,file='plot_energy',status='unknown')
  write(20,*) '# set term x11'
  write(20,*) 'set term postscript landscape monochrome dashed "Helvetica" 22'
  write(20,*)
  write(20,*) 'set xlabel "Time (s)"'
  write(20,*) 'set ylabel "Total energy"'
  write(20,*)
  write(20,*) 'set output "collino3D_total_energy_semilog.eps"'
  write(20,*) 'set logscale y'
  write(20,*) 'plot "energy.dat" t ''Total energy'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)
  close(20)

! create script for Gnuplot
  open(unit=20,file='plotgnu',status='unknown')
  write(20,*) 'set term x11'
  write(20,*) '# set term postscript landscape monochrome dashed "Helvetica" 22'
  write(20,*)
  write(20,*) 'set xlabel "Time (s)"'
  write(20,*) 'set ylabel "Amplitude (m / s)"'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vx_receiver_001.eps"'
  write(20,*) 'plot "Vx_file_001.dat" t ''Vx C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_001.eps"'
  write(20,*) 'plot "Vy_file_001.dat" t ''Vy C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vz_receiver_001.eps"'
  write(20,*) 'plot "Vz_file_001.dat" t ''Vz C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vx_receiver_002.eps"'
  write(20,*) 'plot "Vx_file_002.dat" t ''Vx C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_002.eps"'
  write(20,*) 'plot "Vy_file_002.dat" t ''Vy C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vz_receiver_002.eps"'
  write(20,*) 'plot "Vz_file_002.dat" t ''Vz C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  close(20)

  print *
  print *,'End of the simulation'
  print *

  end program seismic_PML_Collino_3D_iso

!----
!----  save the seismograms in ASCII text format
!----

  subroutine write_seismograms(sisvx,sisvy,nt,nrec,DELTAT)

  implicit none

  integer nt,nrec
  double precision DELTAT

  double precision sisvx(nt,nrec)
  double precision sisvy(nt,nrec)

  integer irec,it

  character(len=100) file_name

! X component
  do irec=1,nrec
    write(file_name,"('Vx_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT),' ',sngl(sisvx(it,irec))
    enddo
    close(11)
  enddo

! Y component
  do irec=1,nrec
    write(file_name,"('Vy_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT),' ',sngl(sisvy(it,irec))
    enddo
    close(11)
  enddo

  end subroutine write_seismograms

!----
!----  routine to create a color image of a given vector component
!----  the image is created in PNM format and then converted to GIF
!----

  subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number)

  implicit none

! non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 0.30d0

! amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0

! use black or white background for points that are below the threshold
  logical, parameter :: WHITE_BACKGROUND = .true.

! size of cross and square in pixels drawn to represent the source and the receivers
  integer, parameter :: width_cross = 5, thickness_cross = 1, size_square = 3

  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX

  double precision, dimension(NX,NY) :: image_data_2D

  integer, dimension(nrec) :: ix_rec,iy_rec

  integer :: ix,iy,irec

  character(len=100) :: file_name,system_command

  integer :: R, G, B

  double precision :: normalized_value,max_amplitude

! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
    write(file_name,"('image',i6.6,'_Vx.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vx.pnm image',i6.6,'_Vx.gif ; rm image',i6.6,'_Vx.pnm')") it,it,it
  else if(field_number == 2) then
    write(file_name,"('image',i6.6,'_Vy.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vy.pnm image',i6.6,'_Vy.gif ; rm image',i6.6,'_Vy.pnm')") it,it,it
  endif

  open(unit=27, file=file_name, status='unknown')

  write(27,"('P3')") ! write image in PNM P3 format

  write(27,*) NX,NY ! write image size
  write(27,*) '255' ! maximum value of each pixel color

! compute maximum amplitude
  max_amplitude = maxval(abs(image_data_2D))

! image starts in upper-left corner in PNM format
  do iy=NY,1,-1
    do ix=1,NX

! define data as vector component normalized to [-1:1] and rounded to nearest integer
! keeping in mind that amplitude can be negative
    normalized_value = image_data_2D(ix,iy) / max_amplitude

! suppress values that are outside [-1:+1] to avoid small edge effects
    if(normalized_value < -1.d0) normalized_value = -1.d0
    if(normalized_value > 1.d0) normalized_value = 1.d0

! draw an orange cross to represent the source
    if((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
        iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
       (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
        iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
      R = 255
      G = 157
      B = 0

! display two-pixel-thick black frame around the image
  else if(ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
      R = 0
      G = 0
      B = 0

! display edges of the PML layers
  else if((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
          (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
          (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
          (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
      R = 255
      G = 150
      B = 0

! suppress all the values that are below the threshold
    else if(abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then

! use a black or white background for points that are below the threshold
      if(WHITE_BACKGROUND) then
        R = 255
        G = 255
        B = 255
      else
        R = 0
        G = 0
        B = 0
      endif

! represent regular image points using red if value is positive, blue if negative
    else if(normalized_value >= 0.d0) then
      R = nint(255.d0*normalized_value**POWER_DISPLAY)
      G = 0
      B = 0
    else
      R = 0
      G = 0
      B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY)
    endif

! draw a green square to represent the receivers
  do irec = 1,nrec
    if((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
        iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square) .or. &
       (ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
        iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square)) then
! use dark green color
      R = 30
      G = 180
      B = 60
    endif
  enddo

! write color pixel
    write(27,"(i3,' ',i3,' ',i3)") R,G,B

    enddo
  enddo

! close file
  close(27)

! call the system to convert image to GIF (can be commented out if "call system" is missing in your compiler)
! call system(system_command)

  end subroutine create_color_image

!
! CeCILL FREE SOFTWARE LICENSE AGREEMENT
!
!     Notice
!
! This Agreement is a Free Software license agreement that is the result
! of discussions between its authors in order to ensure compliance with
! the two main principles guiding its drafting:
!
!     * firstly, compliance with the principles governing the distribution
!       of Free Software: access to source code, broad rights granted to
!       users,
!     * secondly, the election of a governing law, French law, with which
!       it is conformant, both as regards the law of torts and
!       intellectual property law, and the protection that it offers to
!       both authors and holders of the economic rights over software.
!
! The authors of the CeCILL (for Ce[a] C[nrs] I[nria] L[ogiciel] L[ibre])
! license are:
!
! Commissariat a l'Energie Atomique - CEA, a public scientific, technical
! and industrial research establishment, having its principal place of
! business at 25 rue Leblanc, immeuble Le Ponant D, 75015 Paris, France.
!
! Centre National de la Recherche Scientifique - CNRS, a public scientific
! and technological establishment, having its principal place of business
! at 3 rue Michel-Ange, 75794 Paris cedex 16, France.
!
! Institut National de Recherche en Informatique et en Automatique -
! INRIA, a public scientific and technological establishment, having its
! principal place of business at Domaine de Voluceau, Rocquencourt, BP
! 105, 78153 Le Chesnay cedex, France.
!
!     Preamble
!
! The purpose of this Free Software license agreement is to grant users
! the right to modify and redistribute the software governed by this
! license within the framework of an open source distribution model.
!
! The exercising of these rights is conditional upon certain obligations
! for users so as to preserve this status for all subsequent redistributions.
!
! In consideration of access to the source code and the rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors only have limited liability.
!
! In this respect, the risks associated with loading, using, modifying
! and/or developing or reproducing the software by the user are brought to
! the user's attention, given its Free Software status, which may make it
! complicated to use, with the result that its use is reserved for
! developers and experienced professionals having in-depth computer
! knowledge. Users are therefore encouraged to load and test the
! suitability of the software as regards their requirements in conditions
! enabling the security of their systems and/or data to be ensured and,
! more generally, to use and operate it in the same conditions of
! security. This Agreement may be freely reproduced and published,
! provided it is not altered, and that no provisions are either added or
! removed herefrom.
!
! This Agreement may apply to any or all software for which the holder of
! the economic rights decides to submit the use thereof to its provisions.
!
!     Article 1 - DEFINITIONS
!
! For the purpose of this Agreement, when the following expressions
! commence with a capital letter, they shall have the following meaning:
!
! Agreement: means this license agreement, and its possible subsequent
! versions and annexes.
!
! Software: means the software in its Object Code and/or Source Code form
! and, where applicable, its documentation, "as is" when the Licensee
! accepts the Agreement.
!
! Initial Software: means the Software in its Source Code and possibly its
! Object Code form and, where applicable, its documentation, "as is" when
! it is first distributed under the terms and conditions of the Agreement.
!
! Modified Software: means the Software modified by at least one
! Contribution.
!
! Source Code: means all the Software's instructions and program lines to
! which access is required so as to modify the Software.
!
! Object Code: means the binary files originating from the compilation of
! the Source Code.
!
! Holder: means the holder(s) of the economic rights over the Initial
! Software.
!
! Licensee: means the Software user(s) having accepted the Agreement.
!
! Contributor: means a Licensee having made at least one Contribution.
!
! Licensor: means the Holder, or any other individual or legal entity, who
! distributes the Software under the Agreement.
!
! Contribution: means any or all modifications, corrections, translations,
! adaptations and/or new functions integrated into the Software by any or
! all Contributors, as well as any or all Internal Modules.
!
! Module: means a set of sources files including their documentation that
! enables supplementary functions or services in addition to those offered
! by the Software.
!
! External Module: means any or all Modules, not derived from the
! Software, so that this Module and the Software run in separate address
! spaces, with one calling the other when they are run.
!
! Internal Module: means any or all Module, connected to the Software so
! that they both execute in the same address space.
!
! GNU GPL: means the GNU General Public License version 2 or any
! subsequent version, as published by the Free Software Foundation Inc.
!
! Parties: mean both the Licensee and the Licensor.
!
! These expressions may be used both in singular and plural form.
!
!     Article 2 - PURPOSE
!
! The purpose of the Agreement is the grant by the Licensor to the
! Licensee of a non-exclusive, transferable and worldwide license for the
! Software as set forth in Article 5 hereinafter for the whole term of the
! protection granted by the rights over said Software.
!
!     Article 3 - ACCEPTANCE
!
! 3.1 The Licensee shall be deemed as having accepted the terms and
! conditions of this Agreement upon the occurrence of the first of the
! following events:
!
!     * (i) loading the Software by any or all means, notably, by
!       downloading from a remote server, or by loading from a physical
!       medium;
!     * (ii) the first time the Licensee exercises any of the rights
!       granted hereunder.
!
! 3.2 One copy of the Agreement, containing a notice relating to the
! characteristics of the Software, to the limited warranty, and to the
! fact that its use is restricted to experienced users has been provided
! to the Licensee prior to its acceptance as set forth in Article 3.1
! hereinabove, and the Licensee hereby acknowledges that it has read and
! understood it.
!
!     Article 4 - EFFECTIVE DATE AND TERM
!
!       4.1 EFFECTIVE DATE
!
! The Agreement shall become effective on the date when it is accepted by
! the Licensee as set forth in Article 3.1.
!
!       4.2 TERM
!
! The Agreement shall remain in force for the entire legal term of
! protection of the economic rights over the Software.
!
!     Article 5 - SCOPE OF RIGHTS GRANTED
!
! The Licensor hereby grants to the Licensee, who accepts, the following
! rights over the Software for any or all use, and for the term of the
! Agreement, on the basis of the terms and conditions set forth hereinafter.
!
! Besides, if the Licensor owns or comes to own one or more patents
! protecting all or part of the functions of the Software or of its
! components, the Licensor undertakes not to enforce the rights granted by
! these patents against successive Licensees using, exploiting or
! modifying the Software. If these patents are transferred, the Licensor
! undertakes to have the transferees subscribe to the obligations set
! forth in this paragraph.
!
!       5.1 RIGHT OF USE
!
! The Licensee is authorized to use the Software, without any limitation
! as to its fields of application, with it being hereinafter specified
! that this comprises:
!
!    1. permanent or temporary reproduction of all or part of the Software
!       by any or all means and in any or all form.
!
!    2. loading, displaying, running, or storing the Software on any or
!       all medium.
!
!    3. entitlement to observe, study or test its operation so as to
!       determine the ideas and principles behind any or all constituent
!       elements of said Software. This shall apply when the Licensee
!       carries out any or all loading, displaying, running, transmission
!       or storage operation as regards the Software, that it is entitled
!       to carry out hereunder.
!
!       5.2 ENTITLEMENT TO MAKE CONTRIBUTIONS
!
! The right to make Contributions includes the right to translate, adapt,
! arrange, or make any or all modifications to the Software, and the right
! to reproduce the resulting software.
!
! The Licensee is authorized to make any or all Contributions to the
! Software provided that it includes an explicit notice that it is the
! author of said Contribution and indicates the date of the creation thereof.
!
!       5.3 RIGHT OF DISTRIBUTION
!
! In particular, the right of distribution includes the right to publish,
! transmit and communicate the Software to the general public on any or
! all medium, and by any or all means, and the right to market, either in
! consideration of a fee, or free of charge, one or more copies of the
! Software by any means.
!
! The Licensee is further authorized to distribute copies of the modified
! or unmodified Software to third parties according to the terms and
! conditions set forth hereinafter.
!
!         5.3.1 DISTRIBUTION OF SOFTWARE WITHOUT MODIFICATION
!
! The Licensee is authorized to distribute true copies of the Software in
! Source Code or Object Code form, provided that said distribution
! complies with all the provisions of the Agreement and is accompanied by:
!
!    1. a copy of the Agreement,
!
!    2. a notice relating to the limitation of both the Licensor's
!       warranty and liability as set forth in Articles 8 and 9,
!
! and that, in the event that only the Object Code of the Software is
! redistributed, the Licensee allows future Licensees unhindered access to
! the full Source Code of the Software by indicating how to access it, it
! being understood that the additional cost of acquiring the Source Code
! shall not exceed the cost of transferring the data.
!
!         5.3.2 DISTRIBUTION OF MODIFIED SOFTWARE
!
! When the Licensee makes a Contribution to the Software, the terms and
! conditions for the distribution of the resulting Modified Software
! become subject to all the provisions of this Agreement.
!
! The Licensee is authorized to distribute the Modified Software, in
! source code or object code form, provided that said distribution
! complies with all the provisions of the Agreement and is accompanied by:
!
!    1. a copy of the Agreement,
!
!    2. a notice relating to the limitation of both the Licensor's
!       warranty and liability as set forth in Articles 8 and 9,
!
! and that, in the event that only the object code of the Modified
! Software is redistributed, the Licensee allows future Licensees
! unhindered access to the full source code of the Modified Software by
! indicating how to access it, it being understood that the additional
! cost of acquiring the source code shall not exceed the cost of
! transferring the data.
!
!         5.3.3 DISTRIBUTION OF EXTERNAL MODULES
!
! When the Licensee has developed an External Module, the terms and
! conditions of this Agreement do not apply to said External Module, that
! may be distributed under a separate license agreement.
!
!         5.3.4 COMPATIBILITY WITH THE GNU GPL
!
! The Licensee can include a code that is subject to the provisions of one
! of the versions of the GNU GPL in the Modified or unmodified Software,
! and distribute that entire code under the terms of the same version of
! the GNU GPL.
!
! The Licensee can include the Modified or unmodified Software in a code
! that is subject to the provisions of one of the versions of the GNU GPL,
! and distribute that entire code under the terms of the same version of
! the GNU GPL.
!
!     Article 6 - INTELLECTUAL PROPERTY
!
!       6.1 OVER THE INITIAL SOFTWARE
!
! The Holder owns the economic rights over the Initial Software. Any or
! all use of the Initial Software is subject to compliance with the terms
! and conditions under which the Holder has elected to distribute its work
! and no one shall be entitled to modify the terms and conditions for the
! distribution of said Initial Software.
!
! The Holder undertakes that the Initial Software will remain ruled at
! least by this Agreement, for the duration set forth in Article 4.2.
!
!       6.2 OVER THE CONTRIBUTIONS
!
! The Licensee who develops a Contribution is the owner of the
! intellectual property rights over this Contribution as defined by
! applicable law.
!
!       6.3 OVER THE EXTERNAL MODULES
!
! The Licensee who develops an External Module is the owner of the
! intellectual property rights over this External Module as defined by
! applicable law and is free to choose the type of agreement that shall
! govern its distribution.
!
!       6.4 JOINT PROVISIONS
!
! The Licensee expressly undertakes:
!
!    1. not to remove, or modify, in any manner, the intellectual property
!       notices attached to the Software;
!
!    2. to reproduce said notices, in an identical manner, in the copies
!       of the Software modified or not.
!
! The Licensee undertakes not to directly or indirectly infringe the
! intellectual property rights of the Holder and/or Contributors on the
! Software and to take, where applicable, vis-a-vis its staff, any and all
! measures required to ensure respect of said intellectual property rights
! of the Holder and/or Contributors.
!
!     Article 7 - RELATED SERVICES
!
! 7.1 Under no circumstances shall the Agreement oblige the Licensor to
! provide technical assistance or maintenance services for the Software.
!
! However, the Licensor is entitled to offer this type of services. The
! terms and conditions of such technical assistance, and/or such
! maintenance, shall be set forth in a separate instrument. Only the
! Licensor offering said maintenance and/or technical assistance services
! shall incur liability therefor.
!
! 7.2 Similarly, any Licensor is entitled to offer to its licensees, under
! its sole responsibility, a warranty, that shall only be binding upon
! itself, for the redistribution of the Software and/or the Modified
! Software, under terms and conditions that it is free to decide. Said
! warranty, and the financial terms and conditions of its application,
! shall be subject of a separate instrument executed between the Licensor
! and the Licensee.
!
!     Article 8 - LIABILITY
!
! 8.1 Subject to the provisions of Article 8.2, the Licensee shall be
! entitled to claim compensation for any direct loss it may have suffered
! from the Software as a result of a fault on the part of the relevant
! Licensor, subject to providing evidence thereof.
!
! 8.2 The Licensor's liability is limited to the commitments made under
! this Agreement and shall not be incurred as a result of in particular:
! (i) loss due the Licensee's total or partial failure to fulfill its
! obligations, (ii) direct or consequential loss that is suffered by the
! Licensee due to the use or performance of the Software, and (iii) more
! generally, any consequential loss. In particular the Parties expressly
! agree that any or all pecuniary or business loss (i.e. loss of data,
! loss of profits, operating loss, loss of customers or orders,
! opportunity cost, any disturbance to business activities) or any or all
! legal proceedings instituted against the Licensee by a third party,
! shall constitute consequential loss and shall not provide entitlement to
! any or all compensation from the Licensor.
!
!     Article 9 - WARRANTY
!
! 9.1 The Licensee acknowledges that the scientific and technical
! state-of-the-art when the Software was distributed did not enable all
! possible uses to be tested and verified, nor for the presence of
! possible defects to be detected. In this respect, the Licensee's
! attention has been drawn to the risks associated with loading, using,
! modifying and/or developing and reproducing the Software which are
! reserved for experienced users.
!
! The Licensee shall be responsible for verifying, by any or all means,
! the suitability of the product for its requirements, its good working
! order, and for ensuring that it shall not cause damage to either persons
! or properties.
!
! 9.2 The Licensor hereby represents, in good faith, that it is entitled
! to grant all the rights over the Software (including in particular the
! rights set forth in Article 5).
!
! 9.3 The Licensee acknowledges that the Software is supplied "as is" by
! the Licensor without any other express or tacit warranty, other than
! that provided for in Article 9.2 and, in particular, without any warranty
! as to its commercial value, its secured, safe, innovative or relevant
! nature.
!
! Specifically, the Licensor does not warrant that the Software is free
! from any error, that it will operate without interruption, that it will
! be compatible with the Licensee's own equipment and software
! configuration, nor that it will meet the Licensee's requirements.
!
! 9.4 The Licensor does not either expressly or tacitly warrant that the
! Software does not infringe any third party intellectual property right
! relating to a patent, software or any other property right. Therefore,
! the Licensor disclaims any and all liability towards the Licensee
! arising out of any or all proceedings for infringement that may be
! instituted in respect of the use, modification and redistribution of the
! Software. Nevertheless, should such proceedings be instituted against
! the Licensee, the Licensor shall provide it with technical and legal
! assistance for its defense. Such technical and legal assistance shall be
! decided on a case-by-case basis between the relevant Licensor and the
! Licensee pursuant to a memorandum of understanding. The Licensor
! disclaims any and all liability as regards the Licensee's use of the
! name of the Software. No warranty is given as regards the existence of
! prior rights over the name of the Software or as regards the existence
! of a trademark.
!
!     Article 10 - TERMINATION
!
! 10.1 In the event of a breach by the Licensee of its obligations
! hereunder, the Licensor may automatically terminate this Agreement
! thirty (30) days after notice has been sent to the Licensee and has
! remained ineffective.
!
! 10.2 A Licensee whose Agreement is terminated shall no longer be
! authorized to use, modify or distribute the Software. However, any
! licenses that it may have granted prior to termination of the Agreement
! shall remain valid subject to their having been granted in compliance
! with the terms and conditions hereof.
!
!     Article 11 - MISCELLANEOUS
!
!       11.1 EXCUSABLE EVENTS
!
! Neither Party shall be liable for any or all delay, or failure to
! perform the Agreement, that may be attributable to an event of force
! majeure, an act of God or an outside cause, such as defective
! functioning or interruptions of the electricity or telecommunications
! networks, network paralysis following a virus attack, intervention by
! government authorities, natural disasters, water damage, earthquakes,
! fire, explosions, strikes and labor unrest, war, etc.
!
! 11.2 Any failure by either Party, on one or more occasions, to invoke
! one or more of the provisions hereof, shall under no circumstances be
! interpreted as being a waiver by the interested Party of its right to
! invoke said provision(s) subsequently.
!
! 11.3 The Agreement cancels and replaces any or all previous agreements,
! whether written or oral, between the Parties and having the same
! purpose, and constitutes the entirety of the agreement between said
! Parties concerning said purpose. No supplement or modification to the
! terms and conditions hereof shall be effective as between the Parties
! unless it is made in writing and signed by their duly authorized
! representatives.
!
! 11.4 In the event that one or more of the provisions hereof were to
! conflict with a current or future applicable act or legislative text,
! said act or legislative text shall prevail, and the Parties shall make
! the necessary amendments so as to comply with said act or legislative
! text. All other provisions shall remain effective. Similarly, invalidity
! of a provision of the Agreement, for any reason whatsoever, shall not
! cause the Agreement as a whole to be invalid.
!
!       11.5 LANGUAGE
!
! The Agreement is drafted in both French and English and both versions
! are deemed authentic.
!
!     Article 12 - NEW VERSIONS OF THE AGREEMENT
!
! 12.1 Any person is authorized to duplicate and distribute copies of this
! Agreement.
!
! 12.2 So as to ensure coherence, the wording of this Agreement is
! protected and may only be modified by the authors of the License, who
! reserve the right to periodically publish updates or new versions of the
! Agreement, each with a separate number. These subsequent versions may
! address new issues encountered by Free Software.
!
! 12.3 Any Software distributed under a given version of the Agreement may
! only be subsequently distributed under the same version of the Agreement
! or a subsequent version, subject to the provisions of Article 5.3.4.
!
!     Article 13 - GOVERNING LAW AND JURISDICTION
!
! 13.1 The Agreement is governed by French law. The Parties agree to
! endeavor to seek an amicable solution to any disagreements or disputes
! that may arise during the performance of the Agreement.
!
! 13.2 Failing an amicable solution within two (2) months as from their
! occurrence, and unless emergency proceedings are necessary, the
! disagreements or disputes shall be referred to the Paris Courts having
! jurisdiction, by the more diligent Party.
!
! Version 2.0 dated 2006-09-05.
!
