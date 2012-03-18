!
! SEISMIC_CPML Version 1.1.1, November 2009.
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Roland Martin, roland DOT martin aT univ-pau DOT fr
!           and Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the three-dimensional isotropic viscoelastic wave equation
! using a fourth order finite-difference method with Convolutional Perfectly Matched
! Layer (C-PML) conditions.
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

  program seismic_visco_CPML_3D_MPI_OpenMP

! 3D fourth order viscoelastic finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions using 2 mechanisms of attenuation
! with 6 equations per mechanism.

! Roland Martin, University of Pau, France, October 2009.
! based on the elastic code of Komatitsch and Martin, 2007.

! The fourth-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used.

! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
!
! Parallel implementation based on MPI.

! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
! If you use this code for your own research, please cite some (or all) of these
! articles:
!
! @ARTICLE{MaKo09,
! author = {Roland Martin and Dimitri Komatitsch},
! title = {An unsplit convolutional perfectly matched layer technique improved
! at grazing incidence for the viscoelastic wave equation},
! journal = {Geophysical Journal International},
! year = {2009},
! volume = {179},
! pages = {333-344},
! number = {1},
! doi = {10.1111/j.1365-246X.2009.04278.x}}
!
! @ARTICLE{MaKoEz08,
! author = {Roland Martin and Dimitri Komatitsch and Abdela\^aziz Ezziani},
! title = {An unsplit convolutional perfectly matched layer improved at grazing
! incidence for seismic wave equation in poroelastic media},
! journal = {Geophysics},
! year = {2008},
! volume = {73},
! pages = {T51-T61},
! number = {4},
! doi = {10.1190/1.2939484}}
!
! @ARTICLE{MaKoGe08,
! author = {Roland Martin and Dimitri Komatitsch and Stephen D. Gedney},
! title = {A variational formulation of a stabilized unsplit convolutional perfectly
! matched layer for the isotropic or anisotropic seismic wave equation},
! journal = {Computer Modeling in Engineering and Sciences},
! year = {2008},
! volume = {37},
! pages = {274-304},
! number = {3}}
!
! @ARTICLE{KoMa07,
! author = {Dimitri Komatitsch and Roland Martin},
! title = {An unsplit convolutional {P}erfectly {M}atched {L}ayer improved
!          at grazing incidence for the seismic wave equation},
! journal = {Geophysics},
! year = {2007},
! volume = {72},
! number = {5},
! pages = {SM155-SM167},
! doi = {10.1190/1.2757586}}
!
! The original CPML technique for Maxwell's equations is described in:
!
! @ARTICLE{RoGe00,
! author = {J. A. Roden and S. D. Gedney},
! title = {Convolution {PML} ({CPML}): {A}n Efficient {FDTD} Implementation
!          of the {CFS}-{PML} for Arbitrary Media},
! journal = {Microwave and Optical Technology Letters},
! year = {2000},
! volume = {27},
! number = {5},
! pages = {334-339},
! doi = {10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A}}

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

! IMPORTANT : all our CPML codes work fine in single precision as well (which is significantly faster).
!             If you want you can thus force automatic conversion to single precision at compile time
!             or change all the declarations and constants in the code from double precision to single.

  implicit none

! header which contains standard MPI declarations
  include 'mpif.h'

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 210
  integer, parameter :: NY = 800
  integer, parameter :: NZ = 220 ! even number in order to cut along Z axis

! number of processes used in the MPI run
! and local number of points (for simplicity we cut the mesh along Z only)
  integer, parameter :: NPROC = 20
  integer, parameter :: NZ_LOCAL = NZ / NPROC

! size of a grid cell
  double precision, parameter :: DELTAX = 4.d0, ONE_OVER_DELTAX = 1.d0 / DELTAX
  double precision, parameter :: DELTAY = DELTAX, DELTAZ = DELTAX
  double precision, parameter :: ONE_OVER_DELTAY = ONE_OVER_DELTAX, ONE_OVER_DELTAZ = ONE_OVER_DELTAX
  double precision, parameter :: ONE=1.d0,TWO=2.d0, DIM=3.d0
! P-velocity, S-velocity and density
  double precision, parameter :: cp = 3000.d0
  double precision, parameter :: cs = 2000.d0
  double precision, parameter :: rho = 2000.d0
  double precision, parameter :: mu = rho*cs*cs
  double precision, parameter :: lambda = rho*(cp*cp - 2.d0*cs*cs)
  double precision, parameter :: lambdaplustwomu = rho*cp*cp

! total number of time steps
  integer, parameter :: NSTEP = 100000

! time step in seconds
  double precision, parameter :: DELTAT = 4.d-4

! parameters for the source
  double precision, parameter :: f0 = 18.d0
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.d7

! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.
  logical, parameter :: USE_PML_ZMIN = .true.
  logical, parameter :: USE_PML_ZMAX = .true.

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 10

! source
!  integer, parameter :: ISOURCE = NX - 2*NPOINTS_PML - 1
  integer, parameter :: ISOURCE = NPOINTS_PML+20
  integer, parameter :: JSOURCE = NY / 5 + 1
  double precision, parameter :: xsource = (ISOURCE) * DELTAX
  double precision, parameter :: ysource = (JSOURCE) * DELTAY
! angle of source force clockwise with respect to vertical (Y) axis
  double precision, parameter :: ANGLE_FORCE = 0.d0

! receivers
  integer, parameter :: NREC = 3
  double precision, parameter :: xdeb = xsource - 100.d0 ! first receiver x in meters
  double precision, parameter :: ydeb = 2300.d0 ! first receiver y in meters
  double precision, parameter :: xfin = xsource ! last receiver x in meters
  double precision, parameter :: yfin =  300.d0 ! last receiver y in meters

! display information on the screen from time to time
  integer, parameter :: IT_DISPLAY = 10000

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

! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

  double precision, parameter :: K_MAX_PML = 7.d0 ! from Gedney page 8.11
  double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
  double precision, dimension(0:NX+1,0:NY+1,-1:NZ_LOCAL+2) :: &
      memory_dvx_dx, &
      memory_dvx_dy, &
      memory_dvx_dz, &
      memory_dvy_dx, &
      memory_dvy_dy, &
      memory_dvy_dz, &
      memory_dvz_dx, &
      memory_dvz_dy, &
      memory_dvz_dz, &
      memory_dsigmaxx_dx, &
      memory_dsigmayy_dy, &
      memory_dsigmazz_dz, &
      memory_dsigmaxy_dx, &
      memory_dsigmaxy_dy, &
      memory_dsigmaxz_dx, &
      memory_dsigmaxz_dz, &
      memory_dsigmayz_dy, &
      memory_dsigmayz_dz

  double precision :: &
      value_dvx_dx, &
      value_dvx_dy, &
      value_dvx_dz, &
      value_dvy_dx, &
      value_dvy_dy, &
      value_dvy_dz, &
      value_dvz_dx, &
      value_dvz_dy, &
      value_dvz_dz, &
      value_dsigmaxx_dx, &
      value_dsigmayy_dy, &
      value_dsigmazz_dz, &
      value_dsigmaxy_dx, &
      value_dsigmaxy_dy, &
      value_dsigmaxz_dx, &
      value_dsigmaxz_dz, &
      value_dsigmayz_dy, &
      value_dsigmayz_dz

   double precision :: duxdx,duxdy,duxdz,duydx,duydy,duydz,duzdx,duzdy,duzdz,div
! 1D arrays for the damping profiles
  double precision, dimension(1:NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half
  double precision, dimension(1:NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half
  double precision, dimension(1:NZ) :: d_z,K_z,alpha_z,a_z,b_z,d_z_half,K_z_half,alpha_z_half,a_z_half,b_z_half

! PML
  double precision thickness_PML_x,thickness_PML_y,thickness_PML_z
  double precision xoriginleft,xoriginright,yoriginbottom,yorigintop,zoriginbottom,zorigintop
  double precision Rcoef,d0_x,d0_y,d0_z,xval,yval,zval,abscissa_in_PML,abscissa_normalized

! change dimension of Z axis to add two planes for MPI
  double precision, dimension(0:NX+1,0:NY+1,-1:NZ_LOCAL+2) :: vx,vy,vz,sigmaxx,sigmayy,sigmazz,sigmaxy,sigmaxz,sigmayz
  double precision, dimension(0:NX+1,0:NY+1,-1:NZ_LOCAL+2) :: sigmaxx_R,sigmayy_R,sigmazz_R,sigmaxy_R,sigmaxz_R,sigmayz_R
  double precision, dimension(0:NX+1,0:NY+1,-1:NZ_LOCAL+2) :: e1_mech1,e1_mech2,e11_mech1,e11_mech2,e22_mech1,e22_mech2
  double precision, dimension(0:NX+1,0:NY+1,-1:NZ_LOCAL+2) :: e12_mech1,e12_mech2,e13_mech1,e13_mech2,e23_mech1,e23_mech2

  integer, parameter :: number_of_arrays = 9 + 2*9 + 12

! for the source
  double precision a,t,force_x,force_y,source_term

! for receivers
  double precision distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec
  double precision, dimension(NREC) :: xrec,yrec

! for seismograms
  double precision, dimension(NSTEP,NREC) :: sisvx,sisvy

! max amplitude for color snapshots
  double precision max_amplitudeVx
  double precision max_amplitudeVy

! for evolution of total energy in the medium
  double precision :: epsilon_xx,epsilon_yy,epsilon_zz,epsilon_xy,epsilon_xz,epsilon_yz
  double precision, dimension(NSTEP) :: total_energy,total_energy_kinetic,total_energy_potential
  double precision :: local_energy_kinetic,local_energy_potential

  integer :: irec

! precompute some parameters once and for all
  double precision, parameter :: DELTAT_lambda = DELTAT*lambda
  double precision, parameter :: DELTAT_mu = DELTAT*mu
  double precision, parameter :: DELTAT_lambdaplus2mu = DELTAT*lambdaplustwomu

  double precision, parameter :: DELTAT_over_rho = DELTAT/rho
  double precision :: mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed
  double precision :: mul_unrelaxed,lambdal_unrelaxed,lambdalplus2mul_unrelaxed
  double precision :: Un,Sn,Unp1,Mu_nu1,Mu_nu2
  double precision :: phi_nu1_mech1,phi_nu1_mech2
  double precision :: phi_nu2_mech1,phi_nu2_mech2
  double precision :: tauinv,inv_tau_sigma_nu1_mech1,inv_tau_sigma_nu1_mech2
  double precision :: taumin,taumax, tau1, tau2, tau3, tau4
  double precision :: inv_tau_sigma_nu2_mech1,inv_tau_sigma_nu2_mech2
  double precision :: tauinvUn
  double precision :: tau_epsilon_nu1_mech1, tau_sigma_nu1_mech1
  double precision::  tau_epsilon_nu2_mech1, tau_sigma_nu2_mech1
  double precision::  tau_epsilon_nu1_mech2, tau_sigma_nu1_mech2
  double precision::  tau_epsilon_nu2_mech2 ,tau_sigma_nu2_mech2

  integer :: i,j,k,it,it2

  double precision :: Vsolidnorm,Courant_number

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

! array needed for MPI_RECV
  integer, dimension(MPI_STATUS_SIZE) :: message_status

! tag of the message to send
  integer, parameter :: message_tag = 0

! number of values to send or receive
  integer, parameter :: number_of_values = 2*(NX+2)*(NY+2)

  integer :: nb_procs,rank,code,rank_cut_plane,kmin,kmax,kglobal,offset_k,k2begin,kminus1end
  integer :: sender_right_shift,receiver_right_shift,sender_left_shift,receiver_left_shift

!---
!--- program starts here
!---

! start MPI processes
  call MPI_INIT(code)

! get total number of MPI processes in variable nb_procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

! get the rank of our process from 0 (master) to nb_procs-1 (workers)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)

 tau_epsilon_nu1_mech1 = 0.0334d0
  tau_sigma_nu1_mech1   = 0.0303d0

!  tau_epsilon_nu1_mech1 = 0.0325305d0
!  tau_sigma_nu1_mech1   = 0.0311465d0

  tau1= tau_sigma_nu1_mech1/tau_epsilon_nu1_mech1

  tau_epsilon_nu2_mech1 = 0.0352d0
  tau_sigma_nu2_mech1   = 0.0287d0

!  tau_epsilon_nu2_mech1 = 0.0332577d0
!  tau_sigma_nu2_mech1   = 0.0304655d0

  tau2= tau_sigma_nu2_mech1/tau_epsilon_nu2_mech1

  tau_epsilon_nu1_mech2 = 0.0028d0
  tau_sigma_nu1_mech2   = 0.0025d0

!  tau_epsilon_nu1_mech2 = 0.0032530d0
!  tau_sigma_nu1_mech2   = 0.0031146d0

  tau3= tau_sigma_nu1_mech2/tau_epsilon_nu1_mech2

  tau_epsilon_nu2_mech2 = 0.0029d0
  tau_sigma_nu2_mech2   = 0.0024d0

!  tau_epsilon_nu2_mech2 = 0.0033257d0
!  tau_sigma_nu2_mech2   = 0.0030465d0

  tau4= tau_sigma_nu2_mech2/tau_epsilon_nu2_mech2

  taumax=dmax1(1.d0/tau1,1.d0/tau2,1.d0/tau3,1.d0/tau4)
  taumin=dmin1(1.d0/tau1,1.d0/tau2,1.d0/tau3,1.d0/tau4)

 inv_tau_sigma_nu1_mech1 = ONE / tau_sigma_nu1_mech1
  inv_tau_sigma_nu2_mech1 = ONE / tau_sigma_nu2_mech1
  inv_tau_sigma_nu1_mech2 = ONE / tau_sigma_nu1_mech2
  inv_tau_sigma_nu2_mech2 = ONE / tau_sigma_nu2_mech2

phi_nu1_mech1 = (ONE - tau_epsilon_nu1_mech1/tau_sigma_nu1_mech1)&
 / tau_sigma_nu1_mech1
phi_nu2_mech1 = (ONE - tau_epsilon_nu2_mech1/tau_sigma_nu2_mech1)&
 / tau_sigma_nu2_mech1
phi_nu1_mech2 = (ONE - tau_epsilon_nu1_mech2/tau_sigma_nu1_mech2)&
 / tau_sigma_nu1_mech2
phi_nu2_mech2 = (ONE - tau_epsilon_nu2_mech2/tau_sigma_nu2_mech2) &
/ tau_sigma_nu2_mech2

 Mu_nu1 = ONE - (ONE - tau_epsilon_nu1_mech1/tau_sigma_nu1_mech1) &
- (ONE - tau_epsilon_nu1_mech2/tau_sigma_nu1_mech2)
 Mu_nu2 = ONE - (ONE - tau_epsilon_nu2_mech1/tau_sigma_nu2_mech1) &
- (ONE - tau_epsilon_nu2_mech2/tau_sigma_nu2_mech2)

! slice number for the cut plane in the middle of the mesh
  rank_cut_plane = nb_procs/2 - 1

  if(rank == rank_cut_plane) then

  print *
  print *,'3D elastic finite-difference code in velocity and stress formulation with C-PML'
  print *

! display size of the model
  print *
  print *,'NX = ',NX
  print *,'NY = ',NY
  print *,'NZ = ',NZ
  print *
  print *,'NZ_LOCAL = ',NZ_LOCAL
  print *,'NPROC = ',NPROC
  print *
  print *,'size of the model along X = ',(NX+1) * DELTAX
  print *,'size of the model along Y = ',(NY+1) * DELTAY
  print *,'size of the model along Y = ',(NZ+1) * DELTAZ
  print *
  print *,'Total number of grid points = ',(NX+2) * (NY+2) * (NZ+2)
  print *,'Number of points of all the arrays = ',dble(NX+2)*dble(NY+2)*dble(NZ+2)*number_of_arrays
  print *,'Size in GB of all the arrays = ',dble(NX+2)*dble(NY+2)*dble(NZ+2)*number_of_arrays*8.d0/(1024.d0*1024.d0*1024.d0)
  print *
  print *,'In each slice:'
  print *
  print *,'Total number of grid points = ',(NX+2) * (NY+2) * NZ_LOCAL
  print *,'Number of points of the arrays = ',dble(NX+2)*dble(NY+2)*dble(NZ_LOCAL)*number_of_arrays
  print *,'Size in GB of the arrays = ',dble(NX+2)*dble(NY+2)*dble(NZ_LOCAL)*number_of_arrays*8.d0/(1024.d0*1024.d0*1024.d0)
  print *

  endif

! check that code was compiled with the right number of slices
  if(nb_procs /= NPROC) then
    print *,'nb_procs,NPROC = ',nb_procs,NPROC
    stop 'nb_procs must be equal to NPROC'
  endif

! we restrict ourselves to an even number of slices
! in order to have a cut plane in the middle of the mesh for visualization purposes
  if(mod(nb_procs,2) /= 0) stop 'nb_procs must be even'

! check that we can cut along Z in an exact number of slices
  if(mod(NZ,nb_procs) /= 0) stop 'NZ must be a multiple of nb_procs'

! check that a slice is at least as thick as a PML layer
  if(NZ_LOCAL < NPOINTS_PML) stop 'NZ_LOCAL must be greater than NPOINTS_PML'

! offset of this slice when we cut along Z
  offset_k = rank * NZ_LOCAL

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
  thickness_PML_x = NPOINTS_PML * DELTAX
  thickness_PML_y = NPOINTS_PML * DELTAY
  thickness_PML_z = NPOINTS_PML * DELTAZ

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.0001d0

! check that NPOWER is okay
  if(NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * cp *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * cp *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_y)
  d0_z = - (NPOWER + 1) * cp *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_z)

  if(rank == rank_cut_plane) then
    print *
    print *,'d0_x = ',d0_x
    print *,'d0_y = ',d0_y
    print *,'d0_z = ',d0_z
  endif

! PML
  d_x(:) = ZERO
  d_x_half(:) = ZERO
  K_x(:) = 1.d0
  K_x_half(:) = 1.d0
  alpha_x(:) = ZERO
  alpha_x_half(:) = ZERO
  a_x(:) = ZERO
  a_x_half(:) = ZERO

  d_y(:) = ZERO
  d_y_half(:) = ZERO
  K_y(:) = 1.d0
  K_y_half(:) = 1.d0
  alpha_y(:) = ZERO
  alpha_y_half(:) = ZERO
  a_y(:) = ZERO
  a_y_half(:) = ZERO

  d_z(:) = ZERO
  d_z_half(:) = ZERO
  K_z(:) = 1.d0
  K_z_half(:) = 1.d0
  alpha_z(:) = ZERO
  alpha_z_half(:) = ZERO
  a_z(:) = ZERO
  a_z_half(:) = ZERO

! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x
  xoriginright = (NX-1)*DELTAX - thickness_PML_x

  do i = 1,NX

! abscissa of current grid point along the damping profile
    xval = DELTAX * dble(i-1)

!---------- xmin edge
    if(USE_PML_XMIN) then

! define damping profile at the grid points
      abscissa_in_PML = xoriginleft - xval
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0)
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

    endif

!---------- xmax edge
    if(USE_PML_XMAX) then

! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

    endif

! just in case, for -5 at the end
    if(alpha_x(i) < ZERO) alpha_x(i) = ZERO
    if(alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT)

! this to avoid division by zero outside the PML
    if(abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if(abs(d_x_half(i)) > 1.d-6) a_x_half(i) = d_x_half(i) * &
      (b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)))

  enddo

! damping in the Y direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  yoriginbottom = thickness_PML_y
  yorigintop = (NY-1)*DELTAY - thickness_PML_y

  do j = 1,NY

! abscissa of current grid point along the damping profile
    yval = DELTAY * dble(j-1)

!---------- ymin edge
    if(USE_PML_YMIN) then

! define damping profile at the grid points
      abscissa_in_PML = yoriginbottom - yval
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0)
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

    endif

!---------- ymax edge
    if(USE_PML_YMAX) then

! define damping profile at the grid points
      abscissa_in_PML = yval - yorigintop
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

    endif

    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT)
    b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_y_half(j)) * DELTAT)

! this to avoid division by zero outside the PML
    if(abs(d_y(j)) > 1.d-6) a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
    if(abs(d_y_half(j)) > 1.d-6) a_y_half(j) = d_y_half(j) * &
      (b_y_half(j) - 1.d0) / (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_y_half(j)))

  enddo

! damping in the Z direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  zoriginbottom = thickness_PML_z
  zorigintop = (NZ-1)*DELTAZ - thickness_PML_z

  do k = 1,NZ

! abscissa of current grid point along the damping profile
    zval = DELTAZ * dble(k-1)

!---------- zmin edge
    if(USE_PML_ZMIN) then

! define damping profile at the grid points
      abscissa_in_PML = zoriginbottom - zval
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z(k) = d0_z * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

! define damping profile at half the grid points
      abscissa_in_PML = zoriginbottom - (zval + DELTAZ/2.d0)
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z_half(k) = d0_z * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z_half(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z_half(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

    endif

!---------- zmax edge
    if(USE_PML_ZMAX) then

! define damping profile at the grid points
      abscissa_in_PML = zval - zorigintop
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z(k) = d0_z * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

! define damping profile at half the grid points
      abscissa_in_PML = zval + DELTAZ/2.d0 - zorigintop
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z_half(k) = d0_z * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_z_half(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z_half(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) + 0.1d0 * ALPHA_MAX_PML
      endif

    endif

    b_z(k) = exp(- (d_z(k) / K_z(k) + alpha_z(k)) * DELTAT)
    b_z_half(k) = exp(- (d_z_half(k) / K_z_half(k) + alpha_z_half(k)) * DELTAT)

! this to avoid division by zero outside the PML
    if(abs(d_z(k)) > 1.d-6) a_z(k) = d_z(k) * (b_z(k) - 1.d0) / (K_z(k) * (d_z(k) + K_z(k) * alpha_z(k)))
    if(abs(d_z_half(k)) > 1.d-6) a_z_half(k) = d_z_half(k) * &
      (b_z_half(k) - 1.d0) / (K_z_half(k) * (d_z_half(k) + K_z_half(k) * alpha_z_half(k)))

  enddo

  if(rank == rank_cut_plane) then

! print position of the source
  print *
  print *,'Position of the source:'
  print *
  print *,'x = ',xsource
  print *,'y = ',ysource
  print *

! define location of receivers
  print *
  print *,'There are ',nrec,' receivers'
  print *
!  xspacerec = (xfin-xdeb) / dble(NREC-1)
!  yspacerec = (yfin-ydeb) / dble(NREC-1)
!  do irec=1,nrec
!    xrec(irec) = xdeb + dble(irec-1)*xspacerec
!    yrec(irec) = ydeb + dble(irec-1)*yspacerec
!  enddo

 xrec(1)=xsource+500.d0  ! first receiver x in meters
 yrec(1)=ysource+500.d0  ! first receiver y in meters
 xrec(2)=xsource  ! first receiver x in meters
 yrec(2)=ysource+2260.d0  ! first receiver y in meters
 xrec(3)=xsource+500.d0  ! first receiver x in meters
 yrec(3)=ysource+2260.d0  ! first receiver y in meters

! find closest grid point for each receiver
  do irec=1,nrec
    dist = HUGEVAL
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((DELTAX*dble(i) - xrec(irec))**2 + (DELTAY*dble(j) - yrec(irec))**2)
      if(distval < dist) then
        dist = distval
        ix_rec(irec) = i
        iy_rec(irec) = j
      endif
    enddo
    enddo
    print *,'receiver ',irec,' x_target,y_target = ',xrec(irec),yrec(irec)
    print *,'closest grid point found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec)
    print *
  enddo

  endif

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
  Courant_number = cp * dsqrt(taumax)* DELTAT * sqrt(1.d0/DELTAX**2 + 1.d0/DELTAY**2 + 1.d0/DELTAZ**2)
  if(rank == rank_cut_plane) then
    print *,'Courant number is ',Courant_number
    print *,'Vpmax=',cp*dsqrt(taumax)
  endif
  if(Courant_number > 1.d0) stop 'time step is too large, simulation will be unstable'
  print *, "Number of points per wavelength =",cs*dsqrt(taumin)/(2.5d0*f0)/DELTAX,&
   'Vsmin=',cs*dsqrt(taumin)

! erase main arrays
  vx(:,:,:) = ZERO
  vy(:,:,:) = ZERO
  vz(:,:,:) = ZERO

  sigmaxy(:,:,:) = ZERO
  sigmayy(:,:,:) = ZERO
  sigmazz(:,:,:) = ZERO
  sigmaxz(:,:,:) = ZERO
  sigmazz(:,:,:) = ZERO
  sigmayz(:,:,:) = ZERO

  e1_mech1(:,:,:)=ZERO
  e1_mech2(:,:,:)=ZERO
  e11_mech1(:,:,:)=ZERO
  e11_mech2(:,:,:)=ZERO
  e12_mech1(:,:,:)=ZERO
  e12_mech2(:,:,:)=ZERO
  e13_mech1(:,:,:)=ZERO
  e13_mech2(:,:,:)=ZERO
  e23_mech1(:,:,:)=ZERO
  e23_mech2(:,:,:)=ZERO
  e22_mech1(:,:,:)=ZERO
  e22_mech2(:,:,:)=ZERO

! PML
  memory_dvx_dx(:,:,:) = ZERO
  memory_dvx_dy(:,:,:) = ZERO
  memory_dvx_dz(:,:,:) = ZERO
  memory_dvy_dx(:,:,:) = ZERO
  memory_dvy_dy(:,:,:) = ZERO
  memory_dvy_dz(:,:,:) = ZERO
  memory_dvz_dx(:,:,:) = ZERO
  memory_dvz_dy(:,:,:) = ZERO
  memory_dvz_dz(:,:,:) = ZERO
  memory_dsigmaxx_dx(:,:,:) = ZERO
  memory_dsigmayy_dy(:,:,:) = ZERO
  memory_dsigmazz_dz(:,:,:) = ZERO
  memory_dsigmaxy_dx(:,:,:) = ZERO
  memory_dsigmaxy_dy(:,:,:) = ZERO
  memory_dsigmaxz_dx(:,:,:) = ZERO
  memory_dsigmaxz_dz(:,:,:) = ZERO
  memory_dsigmayz_dy(:,:,:) = ZERO
  memory_dsigmayz_dz(:,:,:) = ZERO

! erase seismograms
  sisvx(:,:) = ZERO
  sisvy(:,:) = ZERO

! initialize total energy
  total_energy(:) = ZERO
  total_energy_kinetic(:) = ZERO
  total_energy_potential(:) = ZERO

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

! we receive from the process on the left, and send to the process on the right
  sender_right_shift = rank - 1
  receiver_right_shift = rank + 1

! if we are the first process, there is no neighbor on the left
  if(rank == 0) sender_right_shift = MPI_PROC_NULL

! if we are the last process, there is no neighbor on the right
  if(rank == nb_procs - 1) receiver_right_shift = MPI_PROC_NULL

!---

! we receive from the process on the right, and send to the process on the left
  sender_left_shift = rank + 1
  receiver_left_shift = rank - 1

! if we are the first process, there is no neighbor on the left
  if(rank == 0) receiver_left_shift = MPI_PROC_NULL

! if we are the last process, there is no neighbor on the right
  if(rank == nb_procs - 1) sender_left_shift = MPI_PROC_NULL

  k2begin = 1
  if(rank == 0) k2begin = 2

  kminus1end = NZ_LOCAL
  if(rank == nb_procs - 1) kminus1end = NZ_LOCAL - 1

!---
!---  beginning of time loop
!---

  do it = 1,NSTEP

    if(rank == rank_cut_plane .AND. mod(it,20).eq.0) print *,'it = ',it

!----------------------
! compute stress sigma
!----------------------

! vx(k+1), left shift
  call MPI_SENDRECV(vx(:,:,1:2),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,vx(:,:,NZ_LOCAL+1:NZ_LOCAL+2),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

! vy(k+1), left shift
  call MPI_SENDRECV(vy(:,:,1:2),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,vy(:,:,NZ_LOCAL+1:NZ_LOCAL+2),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

! vz(k-1), right shift
  call MPI_SENDRECV(vz(:,:,NZ_LOCAL-1:NZ_LOCAL),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vz(:,:,-1:0),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  do k=k2begin,NZ_LOCAL
   kglobal = k + offset_k
   do j=2,NY
     do i=1,NX-1

      mul_relaxed = mu
      lambdal_relaxed = lambda
      lambdalplus2mul_relaxed = lambdal_relaxed + TWO*mul_relaxed
      lambdal_unrelaxed = (lambdal_relaxed + 2.d0/DIM*mul_relaxed) * Mu_nu1 - 2.d0/DIM*mul_relaxed * Mu_nu2
      mul_unrelaxed = mul_relaxed * Mu_nu2
      lambdalplus2mul_unrelaxed = lambdal_unrelaxed + TWO*mul_unrelaxed

      value_dvx_dx = (27.d0*vx(i+1,j,k)-27.d0*vx(i,j,k)-vx(i+2,j,k)+vx(i-1,j,k)) * ONE_OVER_DELTAX/24.d0
      value_dvy_dy = (27.d0*vy(i,j,k)-27.d0*vy(i,j-1,k)-vy(i,j+1,k)+vy(i,j-2,k)) * ONE_OVER_DELTAY/24.d0
      value_dvz_dz = (27.d0*vz(i,j,k)-27.d0*vz(i,j,k-1)-vz(i,j,k+1)+vz(i,j,k-2)) * ONE_OVER_DELTAZ/24.d0

      memory_dvx_dx(i,j,k) = b_x_half(i) * memory_dvx_dx(i,j,k) + a_x_half(i) * value_dvx_dx
      memory_dvy_dy(i,j,k) = b_y(j) * memory_dvy_dy(i,j,k) + a_y(j) * value_dvy_dy
      memory_dvz_dz(i,j,k) = b_z(kglobal) * memory_dvz_dz(i,j,k) + a_z(kglobal) * value_dvz_dz

      duxdx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j,k)
      duydy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j,k)
      duzdz = value_dvz_dz / K_z(kglobal) + memory_dvz_dz(i,j,k)

      div=duxdx+duydy+duzdz

!evolution e1_mech1
  tauinv = - inv_tau_sigma_nu1_mech1
  Un = e1_mech1(i,j,k)
  Sn   = div * phi_nu1_mech1
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e1_mech1(i,j,k) = Unp1

!evolution e1_mech2
  tauinv = - inv_tau_sigma_nu1_mech2
  Un = e1_mech2(i,j,k)
  Sn   = div * phi_nu1_mech2
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e1_mech2(i,j,k) = Unp1

! evolution e11_mech1
  tauinv = - inv_tau_sigma_nu2_mech1
  Un = e11_mech1(i,j,k)
  Sn   = (duxdx - div/DIM) * phi_nu2_mech1
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e11_mech1(i,j,k) = Unp1

! evolution e11_mech2
  tauinv = - inv_tau_sigma_nu2_mech2
  Un = e11_mech2(i,j,k)
  Sn   = (duxdx - div/DIM) * phi_nu2_mech2
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e11_mech2(i,j,k) = Unp1

! evolution e22_mech1
  tauinv = - inv_tau_sigma_nu2_mech1
  Un = e22_mech1(i,j,k)
  Sn   = (duydy - div/DIM) * phi_nu2_mech1
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e22_mech1(i,j,k) = Unp1

! evolution e22_mech2
  tauinv = - inv_tau_sigma_nu2_mech2
  Un = e22_mech2(i,j,k)
  Sn   = (duydy - div/DIM) * phi_nu2_mech2
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e22_mech2(i,j,k) = Unp1


!add the memory variables using the relaxed parameters (Carcione page 111)
! : there is a bug in Carcione's equation for sigma_zz
    sigmaxx(i,j,k) = sigmaxx(i,j,k)+deltat*((lambdal_relaxed + 2.d0/DIM*mul_relaxed)* &
      (e1_mech1(i,j,k) + e1_mech2(i,j,k)) + TWO * mul_relaxed * (e11_mech1(i,j,k) + e11_mech2(i,j,k)))
    sigmayy(i,j,k) = sigmayy(i,j,k)+deltat*((lambdal_relaxed + 2.d0/DIM*mul_relaxed)* &
      (e1_mech1(i,j,k) + e1_mech2(i,j,k)) + TWO * mul_relaxed * (e22_mech1(i,j,k) + e22_mech2(i,j,k)))
    sigmazz(i,j,k) = sigmazz(i,j,k)+deltat*((lambdal_relaxed + 2.d0*mul_relaxed)* &
      (e1_mech1(i,j,k) + e1_mech2(i,j,k)) - TWO/DIM * mul_relaxed * (e11_mech1(i,j,k) + e11_mech2(i,j,k)&
      +e22_mech1(i,j,k) + e22_mech2(i,j,k)))

! compute the stress using the unrelaxed Lame parameters (Carcione page 111)

      sigmaxx(i,j,k) = sigmaxx(i,j,k) + &
         (lambdalplus2mul_unrelaxed * (duxdx) + &
          lambdal_unrelaxed* (duydy) + &
          lambdal_unrelaxed* (duzdz) )* DELTAT

      sigmayy(i,j,k) = sigmayy(i,j,k) + &
         (lambdal_unrelaxed * (duxdx) + &
          lambdalplus2mul_unrelaxed* (duydy) +&
          lambdal_unrelaxed* (duzdz)) * DELTAT

      sigmazz(i,j,k) = sigmazz(i,j,k) + &
         (lambdal_unrelaxed * (duxdx) + &
          lambdal_unrelaxed* (duydy) + &
          lambdalplus2mul_unrelaxed* (duzdz)) * DELTAT

      sigmaxx_R(i,j,k) = sigmaxx_R(i,j,k) + &
         (lambdalplus2mul_relaxed * (duxdx) + &
          lambdal_relaxed* (duydy) + &
          lambdal_relaxed* (duzdz) )* DELTAT

      sigmayy_R(i,j,k) = sigmayy_R(i,j,k) + &
         (lambdal_relaxed * (duxdx) + &
          lambdalplus2mul_relaxed* (duydy) +&
          lambdal_relaxed* (duzdz)) * DELTAT

      sigmazz_R(i,j,k) = sigmazz_R(i,j,k) + &
         (lambdal_relaxed * (duxdx) + &
          lambdal_relaxed* (duydy) + &
          lambdalplus2mul_relaxed* (duzdz)) * DELTAT

     enddo
    enddo
  enddo

  do k=1,NZ_LOCAL
   do j=1,NY-1
     do i=2,NX
      mul_relaxed = mu
      mul_unrelaxed = mul_relaxed * Mu_nu2

      value_dvy_dx = (27.d0*vy(i,j,k)-27.d0*vy(i-1,j,k)-vy(i+1,j,k)+vy(i-2,j,k)) * ONE_OVER_DELTAX/24.d0
      value_dvx_dy = (27.d0*vx(i,j+1,k)-27.d0*vx(i,j,k)-vx(i,j+2,k)+vx(i,j-1,k)) * ONE_OVER_DELTAY/24.d0

      memory_dvy_dx(i,j,k) = b_x(i) * memory_dvy_dx(i,j,k) + a_x(i) * value_dvy_dx
      memory_dvx_dy(i,j,k) = b_y_half(j) * memory_dvx_dy(i,j,k) + a_y_half(j) * value_dvx_dy

      duydx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j,k)
      duxdy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j,k)

! evolution e12_mech1
  tauinv = - inv_tau_sigma_nu2_mech1
  Un = e12_mech1(i,j,k)
  Sn   = (duxdy+duydx) * phi_nu2_mech1
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e12_mech1(i,j,k) = Unp1

! evolution e12_mech2
  tauinv = - inv_tau_sigma_nu2_mech2
  Un = e12_mech2(i,j,k)
  Sn   = (duxdy+duydx) * phi_nu2_mech2
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e12_mech2(i,j,k) = Unp1

      sigmaxy(i,j,k) = sigmaxy(i,j,k)+deltat*mul_relaxed * (e12_mech1(i,j,k) + e12_mech2(i,j,k))

    sigmaxy(i,j,k) = sigmaxy(i,j,k) + &
    mul_unrelaxed * (duxdy+duydx) * DELTAT

    sigmaxy_R(i,j,k) = sigmaxy_R(i,j,k) + &
    mul_relaxed * (duxdy+duydx) * DELTAT

      enddo
    enddo
  enddo

  do k=1,kminus1end
   kglobal = k + offset_k
   do j=1,NY
     do i=2,NX
      mul_relaxed = mu
      mul_unrelaxed = mul_relaxed * Mu_nu2

      value_dvz_dx = (27.d0*vz(i,j,k)-27.d0*vz(i-1,j,k)-vz(i+1,j,k)+vz(i-2,j,k)) * ONE_OVER_DELTAX/24.d0
      value_dvx_dz = (27.d0*vx(i,j,k+1)-27.d0*vx(i,j,k)-vx(i,j,k+2)+vx(i,j,k-1)) * ONE_OVER_DELTAZ/24.d0

      memory_dvz_dx(i,j,k) = b_x(i) * memory_dvz_dx(i,j,k) + a_x(i) * value_dvz_dx
      memory_dvx_dz(i,j,k) = b_z_half(kglobal) * memory_dvx_dz(i,j,k) + a_z_half(kglobal) * value_dvx_dz

      duzdx = value_dvz_dx / K_x(i) + memory_dvz_dx(i,j,k)
      duxdz = value_dvx_dz / K_z_half(kglobal) + memory_dvx_dz(i,j,k)

! evolution e13_mech1
  tauinv = - inv_tau_sigma_nu2_mech1
  Un = e13_mech1(i,j,k)
  Sn   = (duxdz+duzdx) * phi_nu2_mech1
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e13_mech1(i,j,k) = Unp1

! evolution e13_mech2
  tauinv = - inv_tau_sigma_nu2_mech2
  Un = e13_mech2(i,j,k)
  Sn   = (duxdz+duzdx) * phi_nu2_mech2
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e13_mech2(i,j,k) = Unp1

      sigmaxz(i,j,k) = sigmaxz(i,j,k)+deltat*mul_relaxed * (e13_mech1(i,j,k) + e13_mech2(i,j,k))

    sigmaxz(i,j,k) = sigmaxz(i,j,k) + &
    mul_unrelaxed * (duxdz+duzdx) * DELTAT

    sigmaxz_R(i,j,k) = sigmaxz_R(i,j,k) + &
    mul_relaxed * (duxdz+duzdx) * DELTAT
      enddo
    enddo

   do j=1,NY-1
     do i=1,NX
      mul_relaxed = mu
      mul_unrelaxed = mul_relaxed * Mu_nu2

      value_dvz_dy = (27.d0*vz(i,j+1,k)-27.d0*vz(i,j,k)-vz(i,j+2,k)+vz(i,j-1,k)) * ONE_OVER_DELTAY/24.d0
      value_dvy_dz = (27.d0*vy(i,j,k+1)-27.d0*vy(i,j,k)-vy(i,j,k+2)+vy(i,j,k-1)) * ONE_OVER_DELTAZ/24.d0

      memory_dvz_dy(i,j,k) = b_y_half(j) * memory_dvz_dy(i,j,k) + a_y_half(j) * value_dvz_dy
      memory_dvy_dz(i,j,k) = b_z_half(kglobal) * memory_dvy_dz(i,j,k) + a_z_half(kglobal) * value_dvy_dz

      duzdy = value_dvz_dy / K_y_half(j) + memory_dvz_dy(i,j,k)
      duydz = value_dvy_dz / K_z_half(kglobal) + memory_dvy_dz(i,j,k)

! evolution e23_mech1
  tauinv = - inv_tau_sigma_nu2_mech1
  Un = e23_mech1(i,j,k)
  Sn   = (duydz+duzdy) * phi_nu2_mech1
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e23_mech1(i,j,k) = Unp1

! evolution e23_mech2
  tauinv = - inv_tau_sigma_nu2_mech2
  Un = e23_mech2(i,j,k)
  Sn   = (duydz+duzdy) * phi_nu2_mech2
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e23_mech2(i,j,k) = Unp1

      sigmayz(i,j,k) = sigmayz(i,j,k)+deltat*mul_relaxed * (e23_mech1(i,j,k) + e23_mech2(i,j,k))

    sigmayz(i,j,k) = sigmayz(i,j,k) + &
    mul_unrelaxed * (duydz+duzdy) * DELTAT

    sigmayz_R(i,j,k) = sigmayz_R(i,j,k) + &
    mul_relaxed * (duydz+duzdy) * DELTAT

      enddo
    enddo
  enddo

!------------------
! compute velocity
!------------------

! sigmazz(k+1), left shift
  call MPI_SENDRECV(sigmazz(:,:,1:2),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,sigmazz(:,:,NZ_LOCAL+1:NZ_LOCAL+2),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

! sigmayz(k-1), right shift
  call MPI_SENDRECV(sigmayz(:,:,NZ_LOCAL-1:NZ_LOCAL),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,sigmayz(:,:,-1:0),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

! sigmaxz(k-1), right shift
  call MPI_SENDRECV(sigmaxz(:,:,NZ_LOCAL-1:NZ_LOCAL),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,sigmaxz(:,:,-1:0),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  do k=k2begin,NZ_LOCAL
   kglobal = k + offset_k
   do j=2,NY
     do i=2,NX

      value_dsigmaxx_dx = (27.d0*sigmaxx(i,j,k)-27.d0*sigmaxx(i-1,j,k)-sigmaxx(i+1,j,k)+sigmaxx(i-2,j,k)) * ONE_OVER_DELTAX/24.d0
      value_dsigmaxy_dy = (27.d0*sigmaxy(i,j,k)-27.d0*sigmaxy(i,j-1,k)-sigmaxy(i,j+1,k)+sigmaxy(i,j-2,k)) * ONE_OVER_DELTAY/24.d0
      value_dsigmaxz_dz = (27.d0*sigmaxz(i,j,k)-27.d0*sigmaxz(i,j,k-1)-sigmaxz(i,j,k+1)+sigmaxz(i,j,k-2)) * ONE_OVER_DELTAZ/24.d0

      memory_dsigmaxx_dx(i,j,k) = b_x(i) * memory_dsigmaxx_dx(i,j,k) + a_x(i) * value_dsigmaxx_dx
      memory_dsigmaxy_dy(i,j,k) = b_y(j) * memory_dsigmaxy_dy(i,j,k) + a_y(j) * value_dsigmaxy_dy
      memory_dsigmaxz_dz(i,j,k) = b_z(kglobal) * memory_dsigmaxz_dz(i,j,k) + a_z(kglobal) * value_dsigmaxz_dz

      value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j,k)
      value_dsigmaxy_dy = value_dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j,k)
      value_dsigmaxz_dz = value_dsigmaxz_dz / K_z(kglobal) + memory_dsigmaxz_dz(i,j,k)

      vx(i,j,k) = DELTAT_over_rho*(value_dsigmaxx_dx + value_dsigmaxy_dy + value_dsigmaxz_dz) + vx(i,j,k)

      enddo
    enddo

   do j=1,NY-1
     do i=1,NX-1

      value_dsigmaxy_dx = (27.d0*sigmaxy(i+1,j,k)-27.d0*sigmaxy(i,j,k)-sigmaxy(i+2,j,k)+sigmaxy(i-1,j,k)) * ONE_OVER_DELTAX/24.d0
      value_dsigmayy_dy = (27.d0*sigmayy(i,j+1,k)-27.d0*sigmayy(i,j,k)-sigmayy(i,j+2,k)+sigmayy(i,j-1,k)) * ONE_OVER_DELTAY/24.d0
      value_dsigmayz_dz = (27.d0*sigmayz(i,j,k)-27.d0*sigmayz(i,j,k-1)-sigmayz(i,j,k+1)+sigmayz(i,j,k-2)) * ONE_OVER_DELTAZ/24.d0

      memory_dsigmaxy_dx(i,j,k) = b_x_half(i) * memory_dsigmaxy_dx(i,j,k) + a_x_half(i) * value_dsigmaxy_dx
      memory_dsigmayy_dy(i,j,k) = b_y_half(j) * memory_dsigmayy_dy(i,j,k) + a_y_half(j) * value_dsigmayy_dy
      memory_dsigmayz_dz(i,j,k) = b_z(kglobal) * memory_dsigmayz_dz(i,j,k) + a_z(kglobal) * value_dsigmayz_dz

      value_dsigmaxy_dx = value_dsigmaxy_dx / K_x_half(i) + memory_dsigmaxy_dx(i,j,k)
      value_dsigmayy_dy = value_dsigmayy_dy / K_y_half(j) + memory_dsigmayy_dy(i,j,k)
      value_dsigmayz_dz = value_dsigmayz_dz / K_z(kglobal) + memory_dsigmayz_dz(i,j,k)

      vy(i,j,k) = DELTAT_over_rho*(value_dsigmaxy_dx + value_dsigmayy_dy + value_dsigmayz_dz) + vy(i,j,k)

      enddo
    enddo
  enddo

  do k=1,kminus1end
   kglobal = k + offset_k
   do j=2,NY
     do i=1,NX-1

      value_dsigmaxz_dx = (27.d0*sigmaxz(i+1,j,k)-27.d0*sigmaxz(i,j,k)-sigmaxz(i+2,j,k)+sigmaxz(i-1,j,k)) * ONE_OVER_DELTAX/24.d0
      value_dsigmayz_dy = (27.d0*sigmayz(i,j,k)-27.d0*sigmayz(i,j-1,k)-sigmayz(i,j+1,k)+sigmayz(i,j-2,k)) * ONE_OVER_DELTAY/24.d0
      value_dsigmazz_dz = (27.d0*sigmazz(i,j,k+1)-27.d0*sigmazz(i,j,k)-sigmazz(i,j,k+2)+sigmazz(i,j,k-1)) * ONE_OVER_DELTAZ/24.d0

      memory_dsigmaxz_dx(i,j,k) = b_x_half(i) * memory_dsigmaxz_dx(i,j,k) + a_x_half(i) * value_dsigmaxz_dx
      memory_dsigmayz_dy(i,j,k) = b_y(j) * memory_dsigmayz_dy(i,j,k) + a_y(j) * value_dsigmayz_dy
      memory_dsigmazz_dz(i,j,k) = b_z_half(kglobal) * memory_dsigmazz_dz(i,j,k) + a_z_half(kglobal) * value_dsigmazz_dz

      value_dsigmaxz_dx = value_dsigmaxz_dx / K_x_half(i) + memory_dsigmaxz_dx(i,j,k)
      value_dsigmayz_dy = value_dsigmayz_dy / K_y(j) + memory_dsigmayz_dy(i,j,k)
      value_dsigmazz_dz = value_dsigmazz_dz / K_z_half(kglobal) + memory_dsigmazz_dz(i,j,k)

      vz(i,j,k) = DELTAT_over_rho*(value_dsigmaxz_dx + value_dsigmayz_dy + value_dsigmazz_dz) + vz(i,j,k)

      enddo
    enddo
  enddo

  if(rank == rank_cut_plane) then

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

! define location of the source
  i = ISOURCE
  j = JSOURCE

  vx(i,j,NZ_LOCAL) = vx(i,j,NZ_LOCAL) + force_x * DELTAT / rho
  vy(i,j,NZ_LOCAL) = vy(i,j,NZ_LOCAL) + force_y * DELTAT / rho

  endif

! implement Dirichlet boundary conditions on the six edges of the grid

! xmin
  vx(0:1,:,:) = ZERO
  vy(0:1,:,:) = ZERO
  vz(0:1,:,:) = ZERO

! xmax
  vx(NX:NX+1,:,:) = ZERO
  vy(NX:NX+1,:,:) = ZERO
  vz(NX:NX+1,:,:) = ZERO

! ymin
  vx(:,0:1,:) = ZERO
  vy(:,0:1,:) = ZERO
  vz(:,0:1,:) = ZERO

! ymax
  vx(:,NY:NY+1,:) = ZERO
  vy(:,NY:NY+1,:) = ZERO
  vz(:,NY:NY+1,:) = ZERO

! zmin
  if(rank == 0) then
    vx(:,:,0:1) = ZERO
    vy(:,:,0:1) = ZERO
    vz(:,:,0:1) = ZERO
  endif

! zmax
  if(rank == nb_procs-1) then
    vx(:,:,NZ_LOCAL:NZ_LOCAL+1) = ZERO
    vy(:,:,NZ_LOCAL:NZ_LOCAL+1) = ZERO
    vz(:,:,NZ_LOCAL:NZ_LOCAL+1) = ZERO
  endif

! store seismograms
  if(rank == rank_cut_plane) then
    do irec = 1,NREC
      sisvx(it,irec) = vx(ix_rec(irec),iy_rec(irec),NZ_LOCAL)
      sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec),NZ_LOCAL)
    enddo
  endif

! compute total energy in the medium (without the PML layers)
  local_energy_kinetic = ZERO
  local_energy_potential = ZERO

  kmin = 1
  kmax = NZ_LOCAL
  if(rank == 0) kmin = NPOINTS_PML
  if(rank == nb_procs-1) kmax = NZ_LOCAL-NPOINTS_PML+1

  do k = kmin,kmax
    do j = NPOINTS_PML, NY-NPOINTS_PML+1
      do i = NPOINTS_PML, NX-NPOINTS_PML+1

! compute kinetic energy first, defined as 1/2 rho ||v||^2
! in principle we should use rho_half_x_half_y instead of rho for vy
! in order to interpolate density at the right location in the staggered grid cell
! but in a homogeneous medium we can safely ignore it
      local_energy_kinetic = local_energy_kinetic + 0.5d0 * rho*( &
              vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2)

! add potential energy, defined as 1/2 epsilon_ij sigma_ij
! in principle we should interpolate the medium parameters at the right location
! in the staggered grid cell but in a homogeneous medium we can safely ignore it

! compute total field from split components
      epsilon_xx = (2.d0*(lambda + mu) * sigmaxx(i,j,k) - lambda * sigmayy(i,j,k) - &
          lambda*sigmazz(i,j,k)) / (2.d0 * mu * (3.d0*lambda + 2.d0*mu))
      epsilon_yy = (2.d0*(lambda + mu) * sigmayy(i,j,k) - lambda * sigmaxx(i,j,k) - &
          lambda*sigmazz(i,j,k)) / (2.d0 * mu * (3.d0*lambda + 2.d0*mu))
      epsilon_zz = (2.d0*(lambda + mu) * sigmazz(i,j,k) - lambda * sigmaxx(i,j,k) - &
          lambda*sigmayy(i,j,k)) / (2.d0 * mu * (3.d0*lambda + 2.d0*mu))
      epsilon_xy = sigmaxy_R(i,j,k) / (2.d0 * mu)
      epsilon_xz = sigmaxz_R(i,j,k) / (2.d0 * mu)
      epsilon_yz = sigmayz_R(i,j,k) / (2.d0 * mu)

      local_energy_potential = local_energy_potential + &
        0.5d0 * (epsilon_xx * sigmaxx_R(i,j,k) + epsilon_yy * sigmayy_R(i,j,k) + &
        epsilon_yy * sigmayy_R(i,j,k)+ 2.d0 * epsilon_xy * sigmaxy_R(i,j,k) + &
        2.d0*epsilon_xz * sigmaxz_R(i,j,k)+2.d0*epsilon_yz * sigmayz_R(i,j,k))

      enddo
    enddo
  enddo

  call MPI_REDUCE(local_energy_kinetic + local_energy_potential,total_energy(it),1, &
                          MPI_DOUBLE_PRECISION,MPI_SUM,rank_cut_plane,MPI_COMM_WORLD,code)
  call MPI_REDUCE(local_energy_kinetic,total_energy_kinetic(it),1, &
                          MPI_DOUBLE_PRECISION,MPI_SUM,rank_cut_plane,MPI_COMM_WORLD,code)
  call MPI_REDUCE(local_energy_potential,total_energy_potential(it),1, &
                          MPI_DOUBLE_PRECISION,MPI_SUM,rank_cut_plane,MPI_COMM_WORLD,code)

! output information
  if(mod(it,IT_DISPLAY) == 0 .or. it == 5) then

    call MPI_REDUCE(maxval(sqrt(vx(:,:,1:NZ_LOCAL)**2 + vy(:,:,1:NZ_LOCAL)**2 + &
        vz(:,:,1:NZ_LOCAL)**2)),Vsolidnorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,rank_cut_plane,MPI_COMM_WORLD,code)

    if(rank == rank_cut_plane) then

      print *,'Time step # ',it
      print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
      print *,'Max norm velocity vector V (m/s) = ',Vsolidnorm
      print *,'Total energy = ',total_energy(it)
! check stability of the code, exit if unstable
      if(Vsolidnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up in solid'

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

! save energy
    open(unit=21,file='energy.dat',status='unknown')
      do it2=1,NSTEP
     write(21,*) sngl(dble(it2-1)*DELTAT),sngl(total_energy_kinetic(it2)),&
        sngl(total_energy_potential(it2)),sngl(total_energy(it2))
      enddo
     close(21)

! save seismograms
    print *,'saving seismograms'
    print *
    call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT,t0)

    call create_color_image(vx(1:NX,1:NY,NZ_LOCAL),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1,max_amplitudeVx)
    call create_color_image(vy(1:NX,1:NY,NZ_LOCAL),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2,max_amplitudeVy)

    endif
    endif

! --- end of time loop
  enddo

  if(rank == rank_cut_plane) then

! save seismograms
  call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT,t0)

! create script for Gnuplot for total energy
  open(unit=20,file='plot_energy',status='unknown')
  write(20,*) '# set term x11'
  write(20,*) 'set term postscript landscape monochrome dashed "Helvetica" 22'
  write(20,*)
  write(20,*) 'set xlabel "Time (s)"'
  write(20,*) 'set ylabel "Total energy"'
  write(20,*)
  write(20,*) 'set output "CPML3D_total_energy_semilog.eps"'
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

  endif

! close MPI program
  call MPI_FINALIZE(code)

  end program seismic_visco_CPML_3D_MPI_OpenMP

!----
!----  save the seismograms in ASCII text format
!----

  subroutine write_seismograms(sisvx,sisvy,nt,nrec,DELTAT,t0)

  implicit none

  integer nt,nrec
  double precision DELTAT,t0

  double precision sisvx(nt,nrec)
  double precision sisvy(nt,nrec)

  integer irec,it

  character(len=100) file_name

! X component
  do irec=1,nrec
    write(file_name,"('Vx_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT-t0),' ',sngl(sisvx(it,irec))
    enddo
    close(11)
  enddo

! Y component
  do irec=1,nrec
    write(file_name,"('Vy_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT-t0),' ',sngl(sisvy(it,irec))
    enddo
    close(11)
  enddo

  end subroutine write_seismograms

!----
!----  routine to create a color image of a given vector component
!----  the image is created in PNM format and then converted to GIF
!----

  subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number,max_amplitude)

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

  character(len=150) :: file_name
! character(len=150) :: system_command

  integer :: R, G, B

  double precision :: normalized_value,max_amplitude

! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org
  if(field_number == 1) then
    write(file_name,"('image',i6.6,'_Vx.pnm')") it
!    write(system_command,"('convert image',i6.6,'_Vx.pnm image',i6.6,'_Vx.gif ; rm image',i6.6,'_Vx.pnm')") it,it,it
  else if(field_number == 2) then
    write(file_name,"('image',i6.6,'_Vy.pnm')") it
!    write(system_command,"('convert image',i6.6,'_Vy.pnm image',i6.6,'_Vy.gif ; rm image',i6.6,'_Vy.pnm')") it,it,it
  endif

  open(unit=27, file=file_name, status='unknown')

  write(27,"('P3')") ! write image in PNM P3 format

  write(27,*) NX,NY ! write image size
  write(27,*) '255' ! maximum value of each pixel color

! compute maximum amplitude
 if(it<=2301) max_amplitude = maxval(abs(image_data_2D))

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
