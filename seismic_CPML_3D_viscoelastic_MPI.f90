!
! SEISMIC_CPML Version 1.2, April 2015.
!
! Copyright CNRS, France.
! Contributors: Roland Martin, roland DOT martin aT get DOT obs-mip DOT fr
!           and Dimitri Komatitsch, komatitsch aT lma DOT cnrs-mrs DOT fr
!
! April 2015: Dimitri Komatitsch added support for the SolvOpt algorithm to compute
! the attenuation parameters in an optimized way. If you use it please cite:
!
! @Article{BlKoChLoXi15,
! Title   = {Positivity-preserving highly-accurate optimization of the {Z}ener viscoelastic model, with application
!            to wave propagation in the presence of strong attenuation},
! Author  = {\'Emilie Blanc and Dimitri Komatitsch and Emmanuel Chaljub and Bruno Lombard and Zhinan Xie},
! Journal = {Geophysical Journal International},
! Year    = {2015},
! Note    = {in press.}}
!
! This software is a computer program whose purpose is to solve
! the three-dimensional isotropic viscoelastic wave equation
! using a fourth order finite-difference method with Convolutional Perfectly Matched Layer (C-PML) conditions.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".

  program seismic_visco_CPML_3D_MPI_OpenMP

! 3D fourth order viscoelastic finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions using 2 mechanisms of attenuation
! with 6 equations per mechanism.

! Roland Martin, University of Pau, France, October 2009.
! based on the elastic code of Komatitsch and Martin, 2007.
! April 2015: Dimitri Komatitsch added support for the SolvOpt algorithm to compute
! the attenuation parameters in an optimized way.

! The fourth-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used.

! *BEWARE* that the attenuation model implemented below is that of J. M. Carcione,
! Seismic modeling in viscoelastic media, Geophysics, vol. 58(1), p. 110-120 (1993), which is NON causal,
! i.e., waves speed up instead of slowing down when turning attenuation on.
! This comes from the fact that in that model the relaxed state at zero frequency is used as a reference instead of
! the unrelaxed state at infinite frequency. These days a causal model should be used instead,
! i.e. one using the unrelaxed state at infinite frequency as a reference.

! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
!
! Parallel implementation based on MPI.

! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
! If you use this code for your own research, please cite some (or all) of these articles:
!
! @Article{BlKoChLoXi15,
! Title   = {Positivity-preserving highly-accurate optimization of the {Z}ener viscoelastic model, with application
!            to wave propagation in the presence of strong attenuation},
! Author  = {\'Emilie Blanc and Dimitri Komatitsch and Emmanuel Chaljub and Bruno Lombard and Zhinan Xie},
! Journal = {Geophysical Journal International},
! Year    = {2015},
! Note    = {in press.}}
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
! doi = {10.1002/1098-2760(20001205)27:5 < 334::AID-MOP14>3.0.CO;2-A}}

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

  use mpi

  implicit none

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 210
  integer, parameter :: NY = 800
  integer, parameter :: NZ = 220 ! even number in order to cut along Z axis

! number of processes used in the MPI run
! and local number of points (for simplicity we cut the mesh along Z only)
  integer, parameter :: NPROC = 4 !! 20
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

! parameters for attenuation
! number of standard linear solids
  integer, parameter :: N_SLS = 2

! Qp approximately equal to 13, Qkappa approximately to 20 and Qmu / Qs approximately to 10
  double precision, parameter :: QKappa_att = 20.d0, QMu_att = 10.d0
  double precision, parameter :: f0_attenuation = 16 ! in Hz

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
! integer, parameter :: ISOURCE = NX - 2*NPOINTS_PML - 1
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

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
  double precision, parameter :: K_MAX_PML = 7.d0
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
  double precision, dimension(N_SLS,0:NX+1,0:NY+1,-1:NZ_LOCAL+2) :: e1,e11,e22,e12,e13,e23

  integer, parameter :: number_of_arrays = 9 + 2*9 + 12

! for the source
  double precision a,t,force_x,force_y,source_term

! for attenuation
  double precision :: f_min_attenuation, f_max_attenuation
  double precision, dimension(N_SLS) :: tau_epsilon_nu1,tau_sigma_nu1,tau_epsilon_nu2,tau_sigma_nu2

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
  double precision :: phi_nu1(N_SLS)
  double precision :: phi_nu2(N_SLS)
  double precision :: tauinv,inv_tau_sigma_nu1(N_SLS)
  double precision :: taumin,taumax,tau1,tau2,tau3,tau4
  double precision :: inv_tau_sigma_nu2(N_SLS)
  double precision :: tauinvUn

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

! attenuation constants for standard linear solids
! nu1 is the dilatation/incompressibility mode (QKappa)
! nu2 is the shear mode (Qmu)
! array index (1) is the first standard linear solid, (2) is the second etc.

! from J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).
! Beware: these values implement specific values of the quality factors:
! Qp approximately equal to 13, Qkappa approximately to 20 and Qmu / Qs approximately to 10,
! which means very high attenuation, see that paper for details.
! tau_epsilon_nu1(1) = 0.0334d0
! tau_sigma_nu1(1)   = 0.0303d0

! tau_epsilon_nu2(1) = 0.0352d0
! tau_sigma_nu2(1)   = 0.0287d0

! tau_epsilon_nu1(2) = 0.0028d0
! tau_sigma_nu1(2)   = 0.0025d0

! tau_epsilon_nu2(2) = 0.0029d0
! tau_sigma_nu2(2)   = 0.0024d0

! from J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).
! Beware: these values implement specific values of the quality factors:
! Qkappa approximately to 27 and Qmu / Qs approximately to 20,
! which means very high attenuation, see that paper for details.

! tau_epsilon_nu1(1) = 0.0325305d0
! tau_sigma_nu1(1)   = 0.0311465d0

! tau_epsilon_nu2(1) = 0.0332577d0
! tau_sigma_nu2(1)   = 0.0304655d0

! tau_epsilon_nu1(2) = 0.0032530d0
! tau_sigma_nu1(2)   = 0.0031146d0

! tau_epsilon_nu2(2) = 0.0033257d0
! tau_sigma_nu2(2)   = 0.0030465d0

! f_min and f_max are computed as : f_max/f_min=12 and (log(f_min)+log(f_max))/2 = log(f0)
  f_min_attenuation = exp(log(f0_attenuation)-log(12.d0)/2.d0)
  f_max_attenuation = 12.d0 * f_min_attenuation

! use new SolvOpt nonlinear optimization with constraints from Emilie Blanc, Bruno Lombard and Dimitri Komatitsch
! to compute attenuation mechanisms
    call compute_attenuation_coeffs(N_SLS,QKappa_att,f0_attenuation,f_min_attenuation,f_max_attenuation, &
                                  tau_epsilon_nu1,tau_sigma_nu1)

    call compute_attenuation_coeffs(N_SLS,QMu_att,f0_attenuation,f_min_attenuation,f_max_attenuation, &
                                  tau_epsilon_nu2,tau_sigma_nu2)

  if (rank == 0) then
    print *
    print *,'with new SolvOpt routine for attenuation:'
    print *
    print *,'N_SLS, QKappa_att, QMu_att = ',N_SLS, QKappa_att, QMu_att
    print *,'f0_attenuation,f_min_attenuation,f_max_attenuation = ',f0_attenuation,f_min_attenuation,f_max_attenuation
    print *,'tau_epsilon_nu1 = ',tau_epsilon_nu1
    print *,'tau_sigma_nu1 = ',tau_sigma_nu1
    print *,'tau_epsilon_nu2 = ',tau_epsilon_nu2
    print *,'tau_sigma_nu2 = ',tau_sigma_nu2
    print *
  endif

  tau1 = tau_sigma_nu1(1)/tau_epsilon_nu1(1)
  tau2 = tau_sigma_nu2(1)/tau_epsilon_nu2(1)
  tau3 = tau_sigma_nu1(2)/tau_epsilon_nu1(2)
  tau4 = tau_sigma_nu2(2)/tau_epsilon_nu2(2)

  taumax = max(1.d0/tau1,1.d0/tau2,1.d0/tau3,1.d0/tau4)
  taumin = min(1.d0/tau1,1.d0/tau2,1.d0/tau3,1.d0/tau4)

  inv_tau_sigma_nu1(1) = ONE / tau_sigma_nu1(1)
  inv_tau_sigma_nu2(1) = ONE / tau_sigma_nu2(1)
  inv_tau_sigma_nu1(2) = ONE / tau_sigma_nu1(2)
  inv_tau_sigma_nu2(2) = ONE / tau_sigma_nu2(2)

  phi_nu1(1) = (ONE - tau_epsilon_nu1(1)/tau_sigma_nu1(1)) / tau_sigma_nu1(1)
  phi_nu2(1) = (ONE - tau_epsilon_nu2(1)/tau_sigma_nu2(1)) / tau_sigma_nu2(1)
  phi_nu1(2) = (ONE - tau_epsilon_nu1(2)/tau_sigma_nu1(2)) / tau_sigma_nu1(2)
  phi_nu2(2) = (ONE - tau_epsilon_nu2(2)/tau_sigma_nu2(2)) / tau_sigma_nu2(2)

  Mu_nu1 = ONE - (ONE - tau_epsilon_nu1(1)/tau_sigma_nu1(1)) - (ONE - tau_epsilon_nu1(2)/tau_sigma_nu1(2))
  Mu_nu2 = ONE - (ONE - tau_epsilon_nu2(1)/tau_sigma_nu2(1)) - (ONE - tau_epsilon_nu2(2)/tau_sigma_nu2(2))

! slice number for the cut plane in the middle of the mesh
  rank_cut_plane = nb_procs/2 - 1

  if (rank == rank_cut_plane) then

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
  if (nb_procs /= NPROC) then
    print *,'error in MPI number of slices: nb_procs,NPROC = ',nb_procs,NPROC,' but they should be equal'
    stop 'nb_procs must be equal to NPROC'
  endif

! we restrict ourselves to an even number of slices
! in order to have a cut plane in the middle of the mesh for visualization purposes
  if (mod(nb_procs,2) /= 0) stop 'nb_procs must be even'

! check that we can cut along Z in an exact number of slices
  if (mod(NZ,nb_procs) /= 0) stop 'NZ must be a multiple of nb_procs'

! check that a slice is at least as thick as a PML layer
  if (NZ_LOCAL < NPOINTS_PML) stop 'NZ_LOCAL must be greater than NPOINTS_PML'

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
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * cp *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * cp *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_y)
  d0_z = - (NPOWER + 1) * cp *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_z)

  if (rank == rank_cut_plane) then
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
    if (USE_PML_XMIN) then

! define damping profile at the grid points
      abscissa_in_PML = xoriginleft - xval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- xmax edge
    if (USE_PML_XMAX) then

! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

! just in case, for -5 at the end
    if (alpha_x(i) < ZERO) alpha_x(i) = ZERO
    if (alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT)

! this to avoid division by zero outside the PML
    if (abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if (abs(d_x_half(i)) > 1.d-6) a_x_half(i) = d_x_half(i) * &
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
    if (USE_PML_YMIN) then

! define damping profile at the grid points
      abscissa_in_PML = yoriginbottom - yval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- ymax edge
    if (USE_PML_YMAX) then

! define damping profile at the grid points
      abscissa_in_PML = yval - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT)
    b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_y_half(j)) * DELTAT)

! this to avoid division by zero outside the PML
    if (abs(d_y(j)) > 1.d-6) a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
    if (abs(d_y_half(j)) > 1.d-6) a_y_half(j) = d_y_half(j) * &
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
    if (USE_PML_ZMIN) then

! define damping profile at the grid points
      abscissa_in_PML = zoriginbottom - zval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z(k) = d0_z * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = zoriginbottom - (zval + DELTAZ/2.d0)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z_half(k) = d0_z * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z_half(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z_half(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- zmax edge
    if (USE_PML_ZMAX) then

! define damping profile at the grid points
      abscissa_in_PML = zval - zorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z(k) = d0_z * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = zval + DELTAZ/2.d0 - zorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_z
        d_z_half(k) = d0_z * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_z_half(k) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_z_half(k) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

    b_z(k) = exp(- (d_z(k) / K_z(k) + alpha_z(k)) * DELTAT)
    b_z_half(k) = exp(- (d_z_half(k) / K_z_half(k) + alpha_z_half(k)) * DELTAT)

! this to avoid division by zero outside the PML
    if (abs(d_z(k)) > 1.d-6) a_z(k) = d_z(k) * (b_z(k) - 1.d0) / (K_z(k) * (d_z(k) + K_z(k) * alpha_z(k)))
    if (abs(d_z_half(k)) > 1.d-6) a_z_half(k) = d_z_half(k) * &
      (b_z_half(k) - 1.d0) / (K_z_half(k) * (d_z_half(k) + K_z_half(k) * alpha_z_half(k)))

  enddo

  if (rank == rank_cut_plane) then

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

! xspacerec = (xfin-xdeb) / dble(NREC-1)
! yspacerec = (yfin-ydeb) / dble(NREC-1)
! do irec=1,nrec
!   xrec(irec) = xdeb + dble(irec-1)*xspacerec
!   yrec(irec) = ydeb + dble(irec-1)*yspacerec
! enddo

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
      if (distval < dist) then
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
  if (rank == rank_cut_plane) then
    print *,'Courant number is ',Courant_number
    print *,'Vpmax=',cp*dsqrt(taumax)
  endif
  if (Courant_number > 1.d0) stop 'time step is too large, simulation will be unstable'
  print *, "Number of points per wavelength =",cs*dsqrt(taumin)/(2.5d0*f0)/DELTAX,'Vsmin=',cs*dsqrt(taumin)

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

  e1(:,:,:,:) = ZERO
  e11(:,:,:,:) = ZERO
  e12(:,:,:,:) = ZERO
  e13(:,:,:,:) = ZERO
  e23(:,:,:,:) = ZERO
  e22(:,:,:,:) = ZERO

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
  if (rank == 0) sender_right_shift = MPI_PROC_NULL

! if we are the last process, there is no neighbor on the right
  if (rank == nb_procs - 1) receiver_right_shift = MPI_PROC_NULL

!---

! we receive from the process on the right, and send to the process on the left
  sender_left_shift = rank + 1
  receiver_left_shift = rank - 1

! if we are the first process, there is no neighbor on the left
  if (rank == 0) receiver_left_shift = MPI_PROC_NULL

! if we are the last process, there is no neighbor on the right
  if (rank == nb_procs - 1) sender_left_shift = MPI_PROC_NULL

  k2begin = 1
  if (rank == 0) k2begin = 2

  kminus1end = NZ_LOCAL
  if (rank == nb_procs - 1) kminus1end = NZ_LOCAL - 1

!---
!---  beginning of time loop
!---

  do it = 1,NSTEP

    if (rank == rank_cut_plane .and. mod(it,20) == 0) print *,'it = ',it

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

! evolution e1(1)
  tauinv = - inv_tau_sigma_nu1(1)
  Un = e1(1,i,j,k)
  Sn   = div * phi_nu1(1)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e1(1,i,j,k) = Unp1

! evolution e1(2)
  tauinv = - inv_tau_sigma_nu1(2)
  Un = e1(2,i,j,k)
  Sn   = div * phi_nu1(2)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e1(2,i,j,k) = Unp1

! evolution e11(1)
  tauinv = - inv_tau_sigma_nu2(1)
  Un = e11(1,i,j,k)
  Sn   = (duxdx - div/DIM) * phi_nu2(1)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e11(1,i,j,k) = Unp1

! evolution e11(2)
  tauinv = - inv_tau_sigma_nu2(2)
  Un = e11(2,i,j,k)
  Sn   = (duxdx - div/DIM) * phi_nu2(2)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e11(2,i,j,k) = Unp1

! evolution e22(1)
  tauinv = - inv_tau_sigma_nu2(1)
  Un = e22(1,i,j,k)
  Sn   = (duydy - div/DIM) * phi_nu2(1)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e22(1,i,j,k) = Unp1

! evolution e22(2)
  tauinv = - inv_tau_sigma_nu2(2)
  Un = e22(2,i,j,k)
  Sn   = (duydy - div/DIM) * phi_nu2(2)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e22(2,i,j,k) = Unp1


!add the memory variables using the relaxed parameters (Carcione page 111)
! : there is a bug in Carcione's equation for sigma_zz
    sigmaxx(i,j,k) = sigmaxx(i,j,k)+deltat*((lambdal_relaxed + 2.d0/DIM*mul_relaxed)* &
      (e1(1,i,j,k) + e1(2,i,j,k)) + TWO * mul_relaxed * (e11(1,i,j,k) + e11(2,i,j,k)))
    sigmayy(i,j,k) = sigmayy(i,j,k)+deltat*((lambdal_relaxed + 2.d0/DIM*mul_relaxed)* &
      (e1(1,i,j,k) + e1(2,i,j,k)) + TWO * mul_relaxed * (e22(1,i,j,k) + e22(2,i,j,k)))
    sigmazz(i,j,k) = sigmazz(i,j,k)+deltat*((lambdal_relaxed + 2.d0*mul_relaxed)* &
      (e1(1,i,j,k) + e1(2,i,j,k)) - TWO/DIM * mul_relaxed * (e11(1,i,j,k) + e11(2,i,j,k)&
      +e22(1,i,j,k) + e22(2,i,j,k)))

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

! evolution e12(1)
  tauinv = - inv_tau_sigma_nu2(1)
  Un = e12(1,i,j,k)
  Sn   = (duxdy+duydx) * phi_nu2(1)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e12(1,i,j,k) = Unp1

! evolution e12(2)
  tauinv = - inv_tau_sigma_nu2(2)
  Un = e12(2,i,j,k)
  Sn   = (duxdy+duydx) * phi_nu2(2)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e12(2,i,j,k) = Unp1

      sigmaxy(i,j,k) = sigmaxy(i,j,k)+deltat*mul_relaxed * (e12(1,i,j,k) + e12(2,i,j,k))

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

! evolution e13(1)
  tauinv = - inv_tau_sigma_nu2(1)
  Un = e13(1,i,j,k)
  Sn   = (duxdz+duzdx) * phi_nu2(1)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e13(1,i,j,k) = Unp1

! evolution e13(2)
  tauinv = - inv_tau_sigma_nu2(2)
  Un = e13(2,i,j,k)
  Sn   = (duxdz+duzdx) * phi_nu2(2)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e13(2,i,j,k) = Unp1

      sigmaxz(i,j,k) = sigmaxz(i,j,k)+deltat*mul_relaxed * (e13(1,i,j,k) + e13(2,i,j,k))

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

! evolution e23(1)
  tauinv = - inv_tau_sigma_nu2(1)
  Un = e23(1,i,j,k)
  Sn   = (duydz+duzdy) * phi_nu2(1)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e23(1,i,j,k) = Unp1

! evolution e23(2)
  tauinv = - inv_tau_sigma_nu2(2)
  Un = e23(2,i,j,k)
  Sn   = (duydz+duzdy) * phi_nu2(2)
  tauinvUn = tauinv * Un
  Unp1 = (Un + deltat*(Sn+0.5d0*tauinvUn))/(1.d0-deltat*0.5d0*tauinv)
  e23(2,i,j,k) = Unp1

      sigmayz(i,j,k) = sigmayz(i,j,k)+deltat*mul_relaxed * (e23(1,i,j,k) + e23(2,i,j,k))

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

  if (rank == rank_cut_plane) then

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
  if (rank == 0) then
    vx(:,:,0:1) = ZERO
    vy(:,:,0:1) = ZERO
    vz(:,:,0:1) = ZERO
  endif

! zmax
  if (rank == nb_procs-1) then
    vx(:,:,NZ_LOCAL:NZ_LOCAL+1) = ZERO
    vy(:,:,NZ_LOCAL:NZ_LOCAL+1) = ZERO
    vz(:,:,NZ_LOCAL:NZ_LOCAL+1) = ZERO
  endif

! store seismograms
  if (rank == rank_cut_plane) then
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
  if (rank == 0) kmin = NPOINTS_PML
  if (rank == nb_procs-1) kmax = NZ_LOCAL-NPOINTS_PML+1

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
  if (mod(it,IT_DISPLAY) == 0 .or. it == 5) then

    call MPI_REDUCE(maxval(sqrt(vx(:,:,1:NZ_LOCAL)**2 + vy(:,:,1:NZ_LOCAL)**2 + &
        vz(:,:,1:NZ_LOCAL)**2)),Vsolidnorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,rank_cut_plane,MPI_COMM_WORLD,code)

    if (rank == rank_cut_plane) then

      print *,'Time step # ',it,' out of ',NSTEP,' out of ',NSTEP
      print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
      print *,'Max norm velocity vector V (m/s) = ',Vsolidnorm
      print *,'Total energy = ',total_energy(it)
! check stability of the code, exit if unstable
      if (Vsolidnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up in solid'

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
     write(21,*) sngl(dble(it2-1)*DELTAT),sngl(total_energy_kinetic(it2)), &
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

  if (rank == rank_cut_plane) then

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

! include the SolvOpt routines
  include "attenuation_model_with_SolvOpt.f90"

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
  if (field_number == 1) then
    write(file_name,"('image',i6.6,'_Vx.pnm')") it
!    write(system_command,"('convert image',i6.6,'_Vx.pnm image',i6.6,'_Vx.gif ; rm image',i6.6,'_Vx.pnm')") it,it,it
  else if (field_number == 2) then
    write(file_name,"('image',i6.6,'_Vy.pnm')") it
!    write(system_command,"('convert image',i6.6,'_Vy.pnm image',i6.6,'_Vy.gif ; rm image',i6.6,'_Vy.pnm')") it,it,it
  endif

  open(unit=27, file=file_name, status='unknown')

  write(27,"('P3')") ! write image in PNM P3 format

  write(27,*) NX,NY ! write image size
  write(27,*) '255' ! maximum value of each pixel color

! compute maximum amplitude
 if (it <= 2301) max_amplitude = maxval(abs(image_data_2D))

! image starts in upper-left corner in PNM format
  do iy=NY,1,-1
    do ix=1,NX

! define data as vector component normalized to [-1:1] and rounded to nearest integer
! keeping in mind that amplitude can be negative
    normalized_value = image_data_2D(ix,iy) / max_amplitude

! suppress values that are outside [-1:+1] to avoid small edge effects
    if (normalized_value < -1.d0) normalized_value = -1.d0
    if (normalized_value > 1.d0) normalized_value = 1.d0

! draw an orange cross to represent the source
    if ((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
        iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
       (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
        iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
      R = 255
      G = 157
      B = 0

! display two-pixel-thick black frame around the image
  else if (ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
      R = 0
      G = 0
      B = 0

! display edges of the PML layers
  else if ((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
          (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
          (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
          (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
      R = 255
      G = 150
      B = 0

! suppress all the values that are below the threshold
    else if (abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then

! use a black or white background for points that are below the threshold
      if (WHITE_BACKGROUND) then
        R = 255
        G = 255
        B = 255
      else
        R = 0
        G = 0
        B = 0
      endif

! represent regular image points using red if value is positive, blue if negative
    else if (normalized_value >= 0.d0) then
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
    if ((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
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

! call the system to convert image to Gif (can be commented out if "call system" is missing in your compiler)
! call system(system_command)

  end subroutine create_color_image

