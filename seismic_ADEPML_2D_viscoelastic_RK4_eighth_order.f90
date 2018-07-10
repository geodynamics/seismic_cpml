!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Roland Martin, roland DOT martin aT get DOT obs-mip DOT fr
!               and Ruiqi Shi and Youshan Liu, China.
!
! RK4 bug detected by Youshan Liu, China fixed by Quentin Brissaud, France and also Caltech (USA) in this version in March 2018.
!
! Ruiqi Shi, Department of Exploration Geophysics, China University of Petroleum, Beijing, China.
! Email: shiruiqi123 AT gmail DOT com
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic wave equation
! using a finite-difference method with Auxiliary Differential
! Equation Perfectly Matched Layer (ADE-PML) conditions.
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

program seismic_ADEPML_2D_viscoelastic_RK4_eighth_order

! High order 2D explicit-semi implicit-implicit viscoelastic finite-difference code
! in velocity and stress formulation with Auxiliary Differential
! Equation Perfectly Matched Layer (ADE-PML) absorbing conditions for
! an SLS viscoelastic medium. It is fourth order Runge-Kutta (RK4) in time
! and 8th order in space using Holberg spatial discretization.

! Version 1.1.3
! by Roland Martin, University of Pau, France, Jan 2010
! with improvements by Ruiqi Shi and
! with a major bug fix in the Runge-Kutta implementation
! and also significant memory usage optimization by Youshan Liu, China, August 2015.
! based on seismic_CPML_2D_isotropic_second_order.f90
! by Dimitri Komatitsch and Roland Martin, University of Pau, France, 2007.

! *BEWARE* that the attenuation model implemented below is that of J. M. Carcione,
! Seismic modeling in viscoelastic media, Geophysics, vol. 58(1), p. 110-120 (1993), which is NON causal,
! i.e., waves speed up instead of slowing down when turning attenuation on.
! This comes from the fact that in that model the relaxed state at zero frequency is used as a reference instead of
! the unrelaxed state at infinite frequency. These days a causal model should be used instead,
! i.e. one using the unrelaxed state at infinite frequency as a reference.

! The 8th-order staggered-grid formulation of Holberg is used:
!
!            ^ y
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |        v_y        |
!   sigma_xy +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!           v_x    sigma_xx
!                  sigma_yy
!

! The ADE-PML implementation is based in part on formulas given in Roden and Gedney (2010)
!
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
! @ARTICLE{MaKoGeBr10,
! author = {Roland Martin and Dimitri Komatitsch and Stephen D. Gedney and Emilien Bruthiaux},
!  title = {A high-order time and space formulation of the unsplit perfectly matched layer
!  for the seismic wave equation using {Auxiliary Differential Equations (ADE-PML)}},
!  journal = {Comput. Model. Eng. Sci.},
!  year = {2010},
!  volume = {56},
!  pages = {17-42},
!  number = {1}}
!
! @ARTICLE{MaCo10,
!  author = {Roland Martin and Carlos Couder-Casta{\~n}eda},
!  title = {An improved unsplit and convolutional Perfectly Matched Layer
!  absorbing technique for the Navier-Stokes equations using cut-off frequency shift},
!  journal = {Comput. Model. Eng. Sci.},
!  pages ={47-77}
!  year = {2010},
!  volume = {63},
!  number = {1}}
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
! @ARTICLE{MaKoEz08,
! author = {Roland Martin and Dimitri Komatitsch and Abdelaaziz Ezziani},
! title = {An unsplit convolutional perfectly matched layer improved at grazing
!          incidence for seismic wave equation in poroelastic media},
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
!          matched layer for the isotropic or anisotropic seismic wave equation},
! journal = {Computer Modeling in Engineering and Sciences},
! year = {2008},
! volume = {37},
! pages = {274-304},
! number = {3}}
!
! @ARTICLE{MaKo09,
!  author = {Roland Martin and Dimitri Komatitsch},
!  title = {An unsplit convolutional perfectly matched layer technique improved
!        at grazing incidence for the viscoelastic wave equation},
!  journal = {Geophysical Journal International},
!  year = {2009},
!  volume = {179},
!  pages = {333-344},
!  number = {1},
!  doi = {10.1111/j.1365-246X.2009.04278.x}}
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
! To display the 2D results as color images, use:
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

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 141
  integer, parameter :: NY = 621  ! NY = 800

! Explicit (epsn=1,epsn=0), implicit (epsn=0,epsn1=1), semi-implicit (epsn=0.5,epsn1=0.5)
  integer, parameter :: iexpl=0
  integer, parameter :: iimpl=0
  integer, parameter :: isemiimpl=1

! size of a grid cell
  double precision, parameter :: DELTAX = 5.d0, ONE_OVER_DELTAX = 1.d0 / DELTAX
  double precision, parameter :: DELTAY = DELTAX
  double precision, parameter :: ONE_OVER_DELTAY = ONE_OVER_DELTAX
  double precision, parameter :: ONE=1.d0,TWO=2.d0, DIM=2.d0

! P-velocity, S-velocity and density
  double precision, parameter :: cp_top = 3050.d0
  double precision, parameter :: cs_top = 1950.d0
  double precision, parameter :: rho_top = 2000.d0
  double precision, parameter :: mu_top = rho_top*cs_top*cs_top
  double precision, parameter :: lambda_top = rho_top*(cp_top*cp_top - 2.d0*cs_top*cs_top)
  double precision, parameter :: lambdaplustwomu_top = rho_top*cp_top*cp_top

  double precision, parameter :: cp_bottom = 2600.d0
  double precision, parameter :: cs_bottom = 1500.d0
  double precision, parameter :: rho_bottom = 1500.d0
  double precision, parameter :: mu_bottom = rho_bottom*cs_bottom*cs_bottom
  double precision, parameter :: lambda_bottom = rho_bottom*(cp_bottom*cp_bottom - 2.d0*cs_bottom*cs_bottom)
  double precision, parameter :: lambdaplustwomu_bottom = rho_bottom*cp_bottom*cp_bottom

! total number of time steps
  integer, parameter :: NSTEP = 5000

! time step in seconds
  double precision, parameter :: DELTAT = 5.d-4

! parameters for the source
  double precision, parameter :: f0 = 15.d0
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.d5

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

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 10

! heterogeneous model and height of the interface
  logical, parameter :: HETEROGENEOUS_MODEL = .true.

! source
! integer, parameter :: ISOURCE = NX - 2*NPOINTS_PML - 1
  integer, parameter :: ISOURCE = NPOINTS_PML+11
  integer, parameter :: JSOURCE = 2*NY / 3
  double precision, parameter :: xsource = (ISOURCE) * DELTAX
  double precision, parameter :: ysource = (JSOURCE) * DELTAY
  double precision, parameter :: INTERFACE_HEIGHT = ysource - 125*DELTAY
  integer, parameter:: JINTERFACE=INT(INTERFACE_HEIGHT/DELTAY)+1
! angle of source force clockwise with respect to vertical (Y) axis
  double precision, parameter :: ANGLE_FORCE = 45.d0

! receivers
  integer, parameter :: NREC = 3
  double precision, parameter :: xdeb = xsource - 100.d0 ! first receiver x in meters
  double precision, parameter :: ydeb = 2300.d0 ! first receiver y in meters
  double precision, parameter :: xfin = xsource ! last receiver x in meters
  double precision, parameter :: yfin =  300.d0 ! last receiver y in meters

! display information on the screen from time to time
  integer, parameter :: IT_DISPLAY = 500

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

! Holberg (1987) coefficients, taken from
! @ARTICLE{Hol87,
! author = {O. Holberg},
! title = {Computational aspects of the choice of operator and sampling interval
! for numerical differentiation in large-scale simulation of wave phenomena},
! journal = {Geophysical Prospecting},
! year = {1987},
! volume = {35},
! pages = {629-655}}
  double precision, parameter :: c1 = 1.231666d0
  double precision, parameter :: c2 = -1.041182d-1
  double precision, parameter :: c3 = 2.063707d-2
  double precision, parameter :: c4 = -3.570998d-3
  double precision, parameter :: coefficient_sum = abs(c1)+abs(c2)+abs(c3)+abs(c4)

! RK4 scheme coefficients, 2 per subloop, 8 in total
  double precision, dimension(4) :: rk41, rk42

! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0
  double precision, parameter :: NPOWER2 = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
  !double precision, parameter :: K_MAX_PML = 7.d0
!  double precision, parameter :: ALPHA_MAX_PML = 0.d0 ! from Festa and Vilotte
  double precision, parameter :: ALPHA_MAX_PML_1 = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte
  double precision K_MAX_PML_1

! double precision, parameter :: K_MAX_PML_2 = K_MAX_PML_1 / 15.d0
!  double precision, parameter :: K_MAX_PML_2 = K_MAX_PML_1
!  double precision, parameter :: ALPHA_MAX_PML_2 =  ALPHA_MAX_PML_1 / 5.d0

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
! We have as many memory variables as the number of frequency shift poles in the CPML
! Indices are 1 and 2 for the 2 frequency shift poles
  double precision, dimension(4,-4:NX+4,-4:NY+4) :: &
      memory_dvx_dx_1, &
      memory_dvx_dy_1, &
      memory_dvy_dx_1, &
      memory_dvy_dy_1, &
      memory_dsigmaxx_dx_1, &
      memory_dsigmayy_dy_1, &
      memory_dsigmaxy_dx_1, &
      memory_dsigmaxy_dy_1
  double precision, dimension(-4:NX+4,-4:NY+4) :: &
      memory_vx_dx_1, &
      memory_vx_dy_1, &
      memory_vy_dx_1, &
      memory_vy_dy_1, &
      memory_sigmaxx_dx_1, &
      memory_sigmayy_dy_1, &
      memory_sigmaxy_dx_1, &
      memory_sigmaxy_dy_1

  double precision :: &
      value_dvx_dx, &
      value_dvx_dy, &
      value_dvy_dx, &
      value_dvy_dy, &
      value_dsigmaxx_dx, &
      value_dsigmayy_dy, &
      value_dsigmaxy_dx, &
      value_dsigmaxy_dy

  double precision :: duxdx,duxdy,duydx,duydy,div
  double precision :: epsn,epsn1,Sn

! 1D arrays for the damping profiles
  double precision, dimension(-4:NX+4) :: d_x_1,K_x_1,alpha_prime_x_1,g_x_1,ksi_x
  double precision, dimension(-4:NX+4) :: d_x_half_1,K_x_half_1,alpha_prime_x_half_1,g_x_half_1,ksi_x_half
  double precision, dimension(-4:NY+4) :: d_y_1,K_y_1,alpha_prime_y_1,g_y_1,ksi_y
  double precision, dimension(-4:NY+4) :: d_y_half_1,K_y_half_1,alpha_prime_y_half_1,g_y_half_1,ksi_y_half

! 1D arrays for the damping profiles
  double precision, dimension(-4:NX+4) :: d_x_2,K_x_2,alpha_prime_x_2,g_x_2
  double precision, dimension(-4:NX+4) :: d_x_half_2,K_x_half_2,alpha_prime_x_half_2,g_x_half_2
  double precision, dimension(-4:NY+4) :: d_y_2,K_y_2,alpha_prime_y_2,g_y_2
  double precision, dimension(-4:NY+4) :: d_y_half_2,K_y_half_2,alpha_prime_y_half_2,g_y_half_2

! coefficients that allow to reset the memory variables at each RK4 substep depend on the substepping and are then of dimension 4,
! 1D arrays for the damping profiles
  double precision, dimension(4,-4:NX+4) :: a_x_1,b_x_1
  double precision, dimension(4,-4:NX+4) :: a_x_half_1,b_x_half_1
  double precision, dimension(4,-4:NY+4) :: a_y_1,b_y_1
  double precision, dimension(4,-4:NY+4) :: a_y_half_1,b_y_half_1

  double precision, dimension(-4:NX+4) :: r_x_1,s_x_1
  double precision, dimension(-4:NX+4) :: r_x_half_1,s_x_half_1
  double precision, dimension(-4:NY+4) :: r_y_1,s_y_1
  double precision, dimension(-4:NY+4) :: r_y_half_1,s_y_half_1

! 1D arrays for the damping profiles
  double precision, dimension(4,-4:NX+4) :: a_x_2
  double precision, dimension(4,-4:NX+4) :: a_x_half_2
  double precision, dimension(4,-4:NY+4) :: a_y_2
  double precision, dimension(4,-4:NY+4) :: a_y_half_2

! PML
  double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
  double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

  double precision, dimension(-4:NX+4,-4:NY+4) :: vx,vy,sigmaxx,sigmayy,sigmaxy
  double precision, dimension(-4:NX+4,-4:NY+4) :: sigmaxx_R,sigmayy_R,sigmaxy_R
  double precision, dimension(N_SLS,-4:NX+4,-4:NY+4) :: e1,e11,e22,e12
  double precision, dimension(-4:NX+4,-4:NY+4) :: rho, mu,lambda,lambdaplustwomu

  double precision rho_half_x_half_y

! variables are stored in four indices in the first dimension to implement RK4
! dv does not always indicate a derivative
  double precision, dimension(4,-4:NX+4,-4:NY+4) :: dvx,dvy,dsigmaxx,dsigmayy,dsigmaxy
  double precision, dimension(4,-4:NX+4,-4:NY+4) :: dsigmaxx_R,dsigmayy_R,dsigmaxy_R
  double precision, dimension(N_SLS,4,-4:NX+4,-4:NY+4) :: de1,de11,de12

  integer, parameter :: number_of_2Darrays = 2*8
  integer, parameter :: number_of_3Darrays = 32

! for the source
  double precision a,t,force_x,force_y,source_term

! for attenuation
  double precision :: f_min_attenuation, f_max_attenuation
  double precision, dimension(N_SLS) :: tau_epsilon_nu1,tau_sigma_nu1,tau_epsilon_nu2,tau_sigma_nu2

! for stability estimate
  double precision :: c_max,c_min

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
  double precision :: epsilon_xx,epsilon_yy,epsilon_xy
  double precision, dimension(NSTEP) :: total_energy,total_energy_kinetic,total_energy_potential
  double precision :: local_energy_kinetic,local_energy_potential

  integer :: irec,inc

  double precision :: mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed
  double precision :: mul_unrelaxed,lambdal_unrelaxed,lambdalplus2mul_unrelaxed
  double precision :: Mu_nu1,Mu_nu2
  double precision :: phi_nu1(N_SLS)
  double precision :: phi_nu2(N_SLS)
  double precision :: tauinv,inv_tau_sigma_nu1(N_SLS)
  double precision :: taumin,taumax, tau1, tau2, tau3, tau4
  double precision :: inv_tau_sigma_nu2(N_SLS)

  integer :: i,j,it,it2

  double precision :: Vsolidnorm

  double precision Courant_number_bottom,Courant_number_top
  double precision Dispersion_number_bottom,Dispersion_number_top

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
!--- the program starts here
!---

  if (iexpl == 1) then
    epsn = 1.d0
    epsn1 = 0.d0
  endif

  if (iimpl == 1) then
    epsn = 0.d0
    epsn1 = 1.d0
  endif

  if (isemiimpl == 1) then
    epsn = 0.5d0
    epsn1 = 0.5d0
  endif

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

  print *
  print *,'2D visco-elastic FD code in velocity and stress formulation with ADE in 8th an RK4'
  print *

! display size of the model
  print *
  print *,'NX = ',NX
  print *,'NY = ',NY
  print *
  print *
  print *,'size of the model along X = ',(NX+1) * DELTAX
  print *,'size of the model along Y = ',(NY+1) * DELTAY
  print *
  print *,'Total number of grid points = ',NX * NY
  print *,'Number of points of all the arrays = ',dble(NX+4*2+1)*dble(NY+4*2+1)*number_of_2Darrays + &
                         4*dble(NX+4*2+1)*dble(NY+4*2+1)*number_of_3Darrays
  print *,'Size in GB of all the arrays = ',dble(NX+4*2+1)*dble(NY+4*2+1)*number_of_2Darrays*8.d0/(1024.d0*1024.d0*1024.d0) + &
                         4*dble(NX+4*2+1)*dble(NY+4*2+1)*number_of_3Darrays*8.d0/(1024.d0*1024.d0*1024.d0)


!--- define profile of absorption in PML region

! thickness of the PML layer in meters
  thickness_PML_x = NPOINTS_PML * DELTAX
  thickness_PML_y = NPOINTS_PML * DELTAY

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
   Rcoef = 1.d-5
  c_max = max(cp_bottom,cp_top)
  c_min = min(cs_bottom,cs_top)

     K_MAX_PML_1 = 1.d0

  print *,'K_MAX_PML = ',K_MAX_PML_1

! check that NPOWER is okay
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  if (HETEROGENEOUS_MODEL) then
  d0_x = - (NPOWER + 1) * c_max *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * c_max *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_y)
 else
  d0_x = - (NPOWER + 1) * cp_bottom *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * cp_bottom *dsqrt(taumax)* log(Rcoef) / (2.d0 * thickness_PML_y)
 endif

    print *
    print *,'d0_x = ',d0_x
    print *,'d0_y = ',d0_y

! parameters involved in RK4 time expansion
  rk41(1) = ZERO
  rk41(2) = 0.5d0
  rk41(3) = 0.5d0
  rk41(4) = 1.d0

  rk42(1) = 1.d0 / 6.d0
  rk42(2) = 2.d0 / 6.d0
  rk42(3) = 2.d0 / 6.d0
  rk42(4) = 1.d0 / 6.d0

  ksi_x(:) = ZERO
  ksi_x_half(:) = ZERO
  d_x_1(:) = ZERO
  d_x_half_1(:) = ZERO
  K_x_1(:) = 1.d0
  K_x_half_1(:) = 1.d0
  alpha_prime_x_1(:) = ZERO
  alpha_prime_x_half_1(:) = ZERO
  a_x_1(:,:) = ZERO
  a_x_half_1(:,:) = ZERO
  g_x_1(:) = 5.d-1
  g_x_half_1(:) = 5.d-1

  ksi_y(:) = ZERO
  ksi_y_half(:) = ZERO
  d_y_1(:) = ZERO
  d_y_half_1(:) = ZERO
  K_y_1(:) = 1.d0
  K_y_half_1(:) = 1.d0
  alpha_prime_y_1(:) = ZERO
  alpha_prime_y_half_1(:) = ZERO
  a_y_1(:,:) = ZERO
  a_y_half_1(:,:) = ZERO
  g_y_1(:) = 1.d0
  g_y_half_1(:) = 1.d0

  d_x_2(:) = ZERO
  d_x_half_2(:) = ZERO
  K_x_2(:) = 1.d0
  K_x_half_2(:) = 1.d0
  alpha_prime_x_2(:) = ZERO
  alpha_prime_x_half_2(:) = ZERO
  a_x_2(:,:) = ZERO
  a_x_half_2(:,:) = ZERO
  g_x_2(:) = 1.d0
  g_x_half_2(:) = 1.d0

  d_y_2(:) = ZERO
  d_y_half_2(:) = ZERO
  K_y_2(:) = 1.d0
  K_y_half_2(:) = 1.d0
  alpha_prime_y_2(:) = ZERO
  alpha_prime_y_half_2(:) = ZERO
  a_y_2(:,:) = ZERO
  a_y_half_2(:,:) = ZERO
  g_y_2(:) = 1.d0
  g_y_half_2(:) =1.d0

  r_x_1(:) = ZERO
  s_x_1(:) = ZERO
  r_x_half_1(:) = ZERO
  s_x_half_1(:) = ZERO
  r_y_1(:) = ZERO
  s_y_1(:) = ZERO
  r_y_half_1(:) = ZERO
  s_y_half_1(:) = ZERO

! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x
  xoriginright = (NX-1)*DELTAX - thickness_PML_x

  do i = -4,NX+4

! abscissa of current grid point along the damping profile
    xval = DELTAX * dble(i-1)

!---------- left edge
    if (USE_PML_XMIN) then

! define damping profile at the grid points
      abscissa_in_PML = xoriginleft - xval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_1(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_1(i) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_x_1(i) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half_1(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half_1(i) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_x_half_1(i) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- right edge
    if (USE_PML_XMAX) then

! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_1(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_1(i) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_x_1(i) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half_1(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half_1(i) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_x_half_1(i) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

    endif

! 1 pole
    d_x_2(i) = 0.d0
    d_x_half_2(i) = 0.d0

! just in case, for -5 at the end
    if (alpha_prime_x_1(i) < ZERO) alpha_prime_x_1(i) = ZERO
    if (alpha_prime_x_half_1(i) < ZERO) alpha_prime_x_half_1(i) = ZERO

! just in case, for -5 at the end
    if (alpha_prime_x_2(i) < ZERO) alpha_prime_x_2(i) = ZERO
    if (alpha_prime_x_half_2(i) < ZERO) alpha_prime_x_half_2(i) = ZERO

! CPML damping parameters for the 4 sub time steps of RK4 algorithm
do inc=1,4
    b_x_1(inc,i) =  (1.-epsn*DELTAT*rk41(inc)*(d_x_1(i)/K_x_1(i) + alpha_prime_x_1(i)))/&
    (1.+epsn1*DELTAT*rk41(inc)*(d_x_1(i)/K_x_1(i) + alpha_prime_x_1(i)))
    b_x_half_1(inc,i) = (1.-epsn*DELTAT*rk41(inc)*(d_x_half_1(i)/K_x_half_1(i) &
   + alpha_prime_x_half_1(i)))/(1. +epsn1*DELTAT*rk41(inc)*(d_x_half_1(i)/K_x_half_1(i) &
    + alpha_prime_x_half_1(i)))

! this to avoid division by zero outside the PML
if (abs(d_x_1(i)) > 1.d-6) a_x_1(inc,i) = - DELTAT*rk41(inc)*d_x_1(i) / (K_x_1(i)* K_x_1(i))/&
 (1. +epsn1*DELTAT*rk41(inc)*(d_x_1(i)/K_x_1(i) + alpha_prime_x_1(i)))

 if (abs(d_x_half_1(i)) > 1.d-6) a_x_half_1(inc,i) =-DELTAT*rk41(inc)*d_x_half_1(i)/&
   (K_x_half_1(i)*K_x_half_1(i) )/&
   (1. +epsn1*DELTAT*rk41(inc)*(d_x_half_1(i)/K_x_half_1(i)&
    + alpha_prime_x_half_1(i)))

   r_x_1(i) = -(d_x_1(i)/K_x_1(i) + alpha_prime_x_1(i))
  s_x_1(i) = - d_x_1(i)/K_x_1(i)/K_x_1(i)
  r_x_half_1(i) = -(d_x_half_1(i)/K_x_half_1(i) + alpha_prime_x_half_1(i))
  s_x_half_1(i) = - d_x_half_1(i)/K_x_half_1(i)/K_x_half_1(i)

  enddo

enddo

! damping in the Y direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  yoriginbottom = thickness_PML_y
  yorigintop = (NY-1)*DELTAY - thickness_PML_y

  do j = -4,NY+4

! abscissa of current grid point along the damping profile
    yval = DELTAY * dble(j-1)

!---------- bottom edge
    if (USE_PML_YMIN) then

! define damping profile at the grid points
      abscissa_in_PML = yoriginbottom - yval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_1(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_1(j) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_y_1(j) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half_1(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half_1(j) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_y_half_1(j) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- top edge
    if (USE_PML_YMAX) then

! define damping profile at the grid points
      abscissa_in_PML = yval - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_1(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_1(j) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_y_1(j) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half_1(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half_1(j) = 1.d0 + (K_MAX_PML_1 - 1.d0) * abscissa_normalized**NPOWER2
        alpha_prime_y_half_1(j) = ALPHA_MAX_PML_1 * (1.d0 - abscissa_normalized)
      endif

    endif

! just in case, for -5 at the end
    if (alpha_prime_y_1(j) < ZERO) alpha_prime_y_1(j) = ZERO
    if (alpha_prime_y_half_1(j) < ZERO) alpha_prime_y_half_1(j) = ZERO

! CPML damping parameters for the 4 sub time steps of RK4 algorithm
do inc=1,4
    b_y_1(inc,j) =  (1.-epsn*DELTAT*rk41(inc)*(d_y_1(j)/K_y_1(j) + alpha_prime_y_1(j)))/&
    (1.+epsn1*DELTAT*rk41(inc)*(d_y_1(j)/K_y_1(j) + alpha_prime_y_1(j)))
    b_y_half_1(inc,j) = (1.-epsn*DELTAT*rk41(inc)*(d_y_half_1(j)/K_y_half_1(j) + &
    alpha_prime_y_half_1(j)))/(1.+epsn1*DELTAT*rk41(inc)*(d_y_half_1(j)/K_y_half_1(j)&
     + alpha_prime_y_half_1(j)))

! this to avoid division by zero outside the PML
  if (abs(d_y_1(j)) > 1.d-6) a_y_1(inc,j) = - DELTAT*rk41(inc)*d_y_1(j)&
   / (K_y_1(j)* K_y_1(j))/&
  (1.+epsn1*DELTAT*rk41(inc)*(d_y_1(j)/K_y_1(j) + alpha_prime_y_1(j)))
 if (abs(d_y_half_1(j)) > 1.d-6) a_y_half_1(inc,j) = -DELTAT*rk41(inc)*d_y_half_1(j) /&
   (K_y_half_1(j) * K_y_half_1(j) )/&
(1.+epsn1*DELTAT*rk41(inc)*(d_y_half_1(j)/K_y_half_1(j) + alpha_prime_y_half_1(j)))
  enddo

  r_y_1(j) = -(d_y_1(j)/K_y_1(j) + alpha_prime_y_1(j))
  s_y_1(j) = - d_y_1(j)/K_y_1(j)/K_y_1(j)
  r_y_half_1(j) = -(d_y_half_1(j)/K_y_half_1(j) + alpha_prime_y_half_1(j))
  s_y_half_1(j) = - d_y_half_1(j)/K_y_half_1(j)/K_y_half_1(j)

enddo

! compute the Lame parameters and density
  do j = -4,NY+4
    do i = -4,NX+4
      if (HETEROGENEOUS_MODEL .and. DELTAY*dble(j-1) > INTERFACE_HEIGHT) then
         rho(i,j)= rho_top
         mu(i,j)= mu_top
         lambda(i,j) = lambda_top
         lambdaplustwomu(i,j) = lambdaplustwomu_top
      else
         rho(i,j)= rho_bottom
         mu(i,j)= mu_bottom
         lambda(i,j) = lambda_bottom
         lambdaplustwomu(i,j) = lambdaplustwomu_bottom
      endif
     enddo
  enddo


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

  xrec(1) = xsource
  yrec(1) = ysource - 393*DELTAY
  xrec(2) = xsource
  yrec(2) = ysource + 191*DELTAY
  xrec(3) = xsource + 101*DELTAX
  yrec(3) = ysource

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

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
  Courant_number_bottom = cp_bottom *dsqrt(taumax)* DELTAT*sqrt(1.d0/DELTAX**2 + 1.d0/DELTAY**2)
  Dispersion_number_bottom=cs_bottom*dsqrt(taumin)/(2.5d0*f0*max(DELTAX,DELTAY))
  print *,'Courant number at the bottom is ',Courant_number_bottom
  print *,'Dispersion number at the bottom is ',Dispersion_number_bottom
  print *
  !if (Courant_number_bottom > 1.d0/coefficient_sum) stop 'time step is too large, simulation will be unstable'

  if (HETEROGENEOUS_MODEL) then
    Courant_number_top = cp_top *dsqrt(taumax) * DELTAT* sqrt(1.d0/DELTAX**2 + 1.d0/DELTAY**2 )
    Dispersion_number_top= cs_top*dsqrt(taumin) /(2.5d0*f0*max(DELTAX,DELTAY))
    print *,'Courant number at the top is ',Courant_number_top
    print *
    print *,'Dispersion number at the top is ',Dispersion_number_top
    !if (Courant_number_top > 1.d0/coefficient_sum) stop 'time step is too large, simulation will be unstable'
  endif

! erase main arrays
  vx(:,:) = ZERO
  vy(:,:) = ZERO
  sigmaxy(:,:) = ZERO
  sigmayy(:,:) = ZERO
  sigmaxx(:,:) = ZERO
  sigmaxy_R(:,:) = ZERO
  sigmayy_R(:,:) = ZERO
  sigmaxx_R(:,:) = ZERO

  dvx(:,:,:) = ZERO
  dvy(:,:,:) = ZERO
  dsigmaxy(:,:,:) = ZERO
  dsigmayy(:,:,:) = ZERO
  dsigmaxx(:,:,:) = ZERO
  dsigmaxy_R(:,:,:) = ZERO
  dsigmayy_R(:,:,:) = ZERO
  dsigmaxx_R(:,:,:) = ZERO

  e1(1,:,:)=ZERO
  e1(2,:,:)=ZERO
  e11(1,:,:)=ZERO
  e11(2,:,:)=ZERO
  e12(1,:,:)=ZERO
  e12(2,:,:)=ZERO
  e22(1,:,:)=ZERO
  e22(2,:,:)=ZERO

  de1(1,:,:,:)=ZERO
  de1(2,:,:,:)=ZERO
  de11(1,:,:,:)=ZERO
  de11(2,:,:,:)=ZERO
  de12(1,:,:,:)=ZERO
  de12(2,:,:,:)=ZERO

! PML
  memory_vx_dx_1(:,:) = ZERO
  memory_vx_dy_1(:,:) = ZERO
  memory_vy_dx_1(:,:) = ZERO
  memory_vy_dy_1(:,:) = ZERO
  memory_sigmaxx_dx_1(:,:) = ZERO
  memory_sigmayy_dy_1(:,:) = ZERO
  memory_sigmaxy_dx_1(:,:) = ZERO
  memory_sigmaxy_dy_1(:,:) = ZERO

  memory_dvx_dx_1(:,:,:) = ZERO
  memory_dvx_dy_1(:,:,:) = ZERO
  memory_dvy_dx_1(:,:,:) = ZERO
  memory_dvy_dy_1(:,:,:) = ZERO
  memory_dsigmaxx_dx_1(:,:,:) = ZERO
  memory_dsigmayy_dy_1(:,:,:) = ZERO
  memory_dsigmaxy_dx_1(:,:,:) = ZERO
  memory_dsigmaxy_dy_1(:,:,:) = ZERO

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
!---  beginning of time loop
!---

  do it = 1,NSTEP
      ! v and sigma temporary variables of RK4
      ! Save initial value for each field

    !dvx(1,:,:) = vx(:,:)
    !dvy(1,:,:) = vy(:,:)
    !dsigmaxx(1,:,:) = sigmaxx(:,:)
    !dsigmayy(1,:,:) = sigmayy(:,:)
    !dsigmaxy(1,:,:) = sigmaxy(:,:)
    !dsigmaxx_R(1,:,:) = sigmaxx_R(:,:)
    !dsigmayy_R(1,:,:) = sigmayy_R(:,:)
    !dsigmaxy_R(1,:,:) = sigmaxy_R(:,:)

    dvx(4,:,:) = vx(:,:)
    dvy(4,:,:) = vy(:,:)
    dsigmaxx(4,:,:) = sigmaxx(:,:)
    dsigmayy(4,:,:) = sigmayy(:,:)
    dsigmaxy(4,:,:) = sigmaxy(:,:)
    dsigmaxx_R(4,:,:) = sigmaxx_R(:,:)
    dsigmayy_R(4,:,:) = sigmayy_R(:,:)
    dsigmaxy_R(4,:,:) = sigmaxy_R(:,:)

    de1(1,4,:,:) = e1(1,:,:)
    de1(2,4,:,:) = e1(2,:,:)
    de11(1,4,:,:) = e11(1,:,:)
    de11(2,4,:,:) = e11(2,:,:)
  ! de22(1,4,:,:) = de22(1,1,:,:)
  ! de22(2,4,:,:) = de22(2,1,:,:)
    de12(1,4,:,:) = e12(1,:,:)
    de12(2,4,:,:) = e12(2,:,:)

    ! same thing for  memory variables
    memory_dsigmaxx_dx_1(4,:,:) = memory_sigmaxx_dx_1(:,:)
    memory_dsigmaxy_dy_1(4,:,:) = memory_sigmaxy_dy_1(:,:)
    memory_dsigmaxy_dx_1(4,:,:) = memory_sigmaxy_dx_1(:,:)
    memory_dsigmayy_dy_1(4,:,:) = memory_sigmayy_dy_1(:,:)
    memory_dvx_dx_1(4,:,:) = memory_vx_dx_1(:,:)
    memory_dvy_dy_1(4,:,:) = memory_vy_dy_1(:,:)
    memory_dvy_dx_1(4,:,:) = memory_vy_dx_1(:,:)
    memory_dvx_dy_1(4,:,:) = memory_vx_dy_1(:,:)

    ! Initialization of time derivatives
    dvx(2,:,:) = ZERO
    dvy(2,:,:) = ZERO
    dsigmaxx(2,:,:) = ZERO
    dsigmayy(2,:,:) = ZERO
    dsigmaxy(2,:,:) = ZERO
    dsigmaxx_R(2,:,:) = ZERO
    dsigmayy_R(2,:,:) = ZERO
    dsigmaxy_R(2,:,:) = ZERO

    de1(1,2,:,:) = ZERO
    de1(2,2,:,:) = ZERO
    de11(1,2,:,:) = ZERO
    de11(2,2,:,:) = ZERO
    de12(1,2,:,:) = ZERO
    de12(2,2,:,:) = ZERO

    ! same thing for  memory variables
    memory_dsigmaxx_dx_1(2,:,:) = ZERO
    memory_dsigmaxy_dy_1(2,:,:) = ZERO
    memory_dsigmaxy_dx_1(2,:,:) = ZERO
    memory_dsigmayy_dy_1(2,:,:) = ZERO
    memory_dvx_dx_1(2,:,:) = ZERO
    memory_dvy_dy_1(2,:,:) = ZERO
    memory_dvy_dx_1(2,:,:) = ZERO
    memory_dvx_dy_1(2,:,:) = ZERO

      ! RK4 loop (loop on the four RK4 substeps)
    do inc= 1,4

! The new values of the different variables v and sigma are computed
        dvx(1,:,:) = dvx(4,:,:) + rk41(inc) * dvx(2,:,:) * DELTAT
        dvy(1,:,:) = dvy(4,:,:) + rk41(inc) * dvy(2,:,:) * DELTAT
        dsigmaxx(1,:,:) = dsigmaxx(4,:,:) + rk41(inc) * dsigmaxx(2,:,:) * DELTAT
        dsigmayy(1,:,:) = dsigmayy(4,:,:) + rk41(inc) * dsigmayy(2,:,:) * DELTAT
        dsigmaxy(1,:,:) = dsigmaxy(4,:,:) + rk41(inc) * dsigmaxy(2,:,:) * DELTAT
        dsigmaxx_R(1,:,:) = dsigmaxx_R(4,:,:) + rk41(inc) * dsigmaxx_R(2,:,:) * DELTAT
        dsigmayy_R(1,:,:) = dsigmayy_R(4,:,:) + rk41(inc) * dsigmayy_R(2,:,:) * DELTAT
        dsigmaxy_R(1,:,:) = dsigmaxy_R(4,:,:) + rk41(inc) * dsigmaxy_R(2,:,:) * DELTAT

        de1(1,1,:,:) = de1(1,4,:,:) + rk41(inc) * de1(1,2,:,:) * DELTAT
        de1(2,1,:,:) = de1(2,4,:,:) + rk41(inc) * de1(2,2,:,:) * DELTAT
        de11(1,1,:,:) = de11(1,4,:,:) + rk41(inc) * de11(1,2,:,:) * DELTAT
        de11(2,1,:,:) = de11(2,4,:,:) + rk41(inc) * de11(2,2,:,:) * DELTAT
        de12(1,1,:,:) = de12(1,4,:,:) + rk41(inc) * de12(1,2,:,:) * DELTAT
        de12(2,1,:,:) = de12(2,4,:,:) + rk41(inc) * de12(2,2,:,:) * DELTAT

        memory_dsigmaxx_dx_1(1,:,:) = memory_dsigmaxx_dx_1(4,:,:) + rk41(inc)*DELTAT*memory_dsigmaxx_dx_1(2,:,:)
        memory_dsigmaxy_dy_1(1,:,:) = memory_dsigmaxy_dy_1(4,:,:) + rk41(inc)*DELTAT*memory_dsigmaxy_dy_1(2,:,:)
        memory_dsigmaxy_dx_1(1,:,:) = memory_dsigmaxy_dx_1(4,:,:) + rk41(inc)*DELTAT*memory_dsigmaxy_dx_1(2,:,:)
        memory_dsigmayy_dy_1(1,:,:) = memory_dsigmayy_dy_1(4,:,:) + rk41(inc)*DELTAT*memory_dsigmayy_dy_1(2,:,:)
        memory_dvx_dx_1(1,:,:) = memory_dvx_dx_1(4,:,:) + rk41(inc)*DELTAT*memory_dvx_dx_1(2,:,:)
        memory_dvy_dy_1(1,:,:) = memory_dvy_dy_1(4,:,:) + rk41(inc)*DELTAT*memory_dvy_dy_1(2,:,:)
        memory_dvx_dy_1(1,:,:) = memory_dvx_dy_1(4,:,:) + rk41(inc)*DELTAT*memory_dvx_dy_1(2,:,:)
        memory_dvy_dx_1(1,:,:) = memory_dvy_dx_1(4,:,:) + rk41(inc)*DELTAT*memory_dvy_dx_1(2,:,:)

     !------------------
     ! compute velocity
     !------------------
      do j = 2,NY
            do i = 2,NX

          value_dsigmaxx_dx = ( c1 * (dsigmaxx(1,i,j) - dsigmaxx(1,i-1,j)) + c2 * (dsigmaxx(1,i+1,j) - dsigmaxx(1,i-2,j)) + &
                    c3 * (dsigmaxx(1,i+2,j) - dsigmaxx(1,i-3,j)) + c4 * (dsigmaxx(1,i+3,j) - dsigmaxx(1,i-4,j)) ) * ONE_OVER_DELTAX

          value_dsigmaxy_dy = ( c1 * (dsigmaxy(1,i,j) - dsigmaxy(1,i,j-1)) + c2* (dsigmaxy(1,i,j+1) - dsigmaxy(1,i,j-2)) + &
                    c3 * (dsigmaxy(1,i,j+2) - dsigmaxy(1,i,j-3)) + c4 * (dsigmaxy(1,i,j+3) - dsigmaxy(1,i,j-4)) ) * ONE_OVER_DELTAY

          if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then

          memory_dsigmaxx_dx_1(2,i,j) = r_x_1(i) * memory_dsigmaxx_dx_1(1,i,j) + s_x_1(i) * value_dsigmaxx_dx
          memory_dsigmaxy_dy_1(2,i,j) = r_y_1(j) * memory_dsigmaxy_dy_1(1,i,j) + s_y_1(j) * value_dsigmaxy_dy

        value_dsigmaxx_dx = value_dsigmaxx_dx / K_x_1(i) + memory_dsigmaxx_dx_1(1,i,j)
          value_dsigmaxy_dy = value_dsigmaxy_dy / K_y_1(j) + memory_dsigmaxy_dy_1(1,i,j)
          endif

          dvx(2,i,j) = (value_dsigmaxx_dx + value_dsigmaxy_dy)/rho(i,j)

            enddo
        enddo

        do j = 1,NY-1
            do i = 1,NX-1
             rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))

             value_dsigmaxy_dx = ( c1 * (dsigmaxy(1,i+1,j) - dsigmaxy(1,i,j)) + c2 * (dsigmaxy(1,i+2,j) - dsigmaxy(1,i-1,j)) + &
                    c3 * (dsigmaxy(1,i+3,j) - dsigmaxy(1,i-2,j)) + c4 * (dsigmaxy(1,i+4,j) - dsigmaxy(1,i-3,j)) )* ONE_OVER_DELTAX

             value_dsigmayy_dy = ( c1 * (dsigmayy(1,i,j+1) - dsigmayy(1,i,j)) + c2 * (dsigmayy(1,i,j+2) - dsigmayy(1,i,j-1)) + &
                    c3 * (dsigmayy(1,i,j+3) - dsigmayy(1,i,j-2)) + c4 * (dsigmayy(1,i,j+4) - dsigmayy(1,i,j-3)) )* ONE_OVER_DELTAY

            if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then
            memory_dsigmaxy_dx_1(2,i,j) = r_x_half_1(i) * memory_dsigmaxy_dx_1(1,i,j) + s_x_half_1(i) * value_dsigmaxy_dx
            memory_dsigmayy_dy_1(2,i,j) = r_y_half_1(j) * memory_dsigmayy_dy_1(1,i,j) + s_y_half_1(j) * value_dsigmayy_dy

            value_dsigmaxy_dx = value_dsigmaxy_dx/K_x_half_1(i)+memory_dsigmaxy_dx_1(1,i,j)
            value_dsigmayy_dy = value_dsigmayy_dy/K_y_half_1(j)+memory_dsigmayy_dy_1(1,i,j)
            endif

                dvy(2,i,j) = (value_dsigmaxy_dx + value_dsigmayy_dy) /rho_half_x_half_y
            enddo
        enddo


    ! add the source (force vector located at a given grid point)
     a = pi*pi*f0*f0;
     t = (dble(it-1)+ rk41(inc)) * DELTAT

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

    ! interpolate density at the right location in the staggered grid cell
    dvx(2,i,j) = dvx(2,i,j) + force_x/ rho(i,j)

    rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))
    dvy(2,i,j) = dvy(2,i,j) + force_y/ rho_half_x_half_y

    ! Dirichlet conditions (rigid boundaries) on all the edges of the grid
        dvx(:,-4:1,:) = ZERO
        dvx(:,NX:NX+4,:) = ZERO

        dvx(:,:,-4:1) = ZERO
        dvx(:,:,NY:NY+4) = ZERO

        dvy(:,-4:1,:) = ZERO
        dvy(:,NX:NX+4,:) = ZERO

        dvy(:,:,-4:1) = ZERO
        dvy(:,:,NY:NY+4) = ZERO

   !----------------------
   ! compute stress sigma
   !----------------------

   do j=2,NY
     do i=1,NX-1

      mul_relaxed = 0.5d0 * (mu(i+1,j) + mu(i,j))
      lambdal_relaxed = 0.5d0 * (lambda(i+1,j) + lambda(i,j))
      lambdalplus2mul_relaxed = 0.5d0 * (lambdaplustwomu(i+1,j) + lambdaplustwomu(i,j))

      lambdal_unrelaxed = (lambdal_relaxed + 2.d0/DIM*mul_relaxed) * Mu_nu1 - 2.d0/DIM*mul_relaxed * Mu_nu2
      mul_unrelaxed = mul_relaxed * Mu_nu2
      lambdalplus2mul_unrelaxed = lambdal_unrelaxed + TWO*mul_unrelaxed

        value_dvx_dx = ( c1 * (dvx(1,i+1,j) - dvx(1,i,j)) + c2 * (dvx(1,i+2,j) - dvx(1,i-1,j)) + &
                  c3 * (dvx(1,i+3,j) - dvx(1,i-2,j)) + c4 * (dvx(1,i+4,j) - dvx(1,i-3,j)) )* ONE_OVER_DELTAX

        value_dvy_dy = ( c1 * (dvy(1,i,j) - dvy(1,i,j-1)) + c2 * (dvy(1,i,j+1) - dvy(1,i,j-2)) + &
                  c3 * (dvy(1,i,j+2) - dvy(1,i,j-3)) + c4 * (dvy(1,i,j+3) - dvy(1,i,j-4)) )* ONE_OVER_DELTAY

        duxdx = value_dvx_dx
        duydy = value_dvy_dy

      if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then
           memory_dvx_dx_1(2,i,j) = r_x_half_1(i) * memory_dvx_dx_1(1,i,j) + s_x_half_1(i) * value_dvx_dx
           memory_dvy_dy_1(2,i,j) = r_y_1(j) * memory_dvy_dy_1(1,i,j) + s_y_1(j) * value_dvy_dy

           duxdx = value_dvx_dx / K_x_half_1(i) + memory_dvx_dx_1(1,i,j)
           duydy = value_dvy_dy / K_y_1(j) + memory_dvy_dy_1(1,i,j)
        endif

      div=duxdx+duydy

!evolution e1(1)
 tauinv = - inv_tau_sigma_nu1(1)
 Sn   = div * phi_nu1(1)
 de1(1,2,i,j) = tauinv * de1(1,1,i,j) + Sn

!evolution e1(2)
 tauinv = - inv_tau_sigma_nu1(2)
 Sn   = div * phi_nu1(2)
 de1(2,2,i,j) = tauinv * de1(2,1,i,j) + Sn

! evolution e11(1)
 tauinv = - inv_tau_sigma_nu2(1)
 Sn   = (duxdx - div/DIM) * phi_nu2(1)
 de11(1,2,i,j) = tauinv * de11(1,1,i,j) + Sn

! evolution e11(2)
 tauinv = - inv_tau_sigma_nu2(2)
 Sn   = (duxdx - div/DIM) * phi_nu2(2)
 de11(2,2,i,j) = tauinv * de11(2,1,i,j) + Sn

!add the memory variables using the relaxed parameters (Carcione page 111)
! there is a bug in Carcione's equation for sigma_zz
  dsigmaxx(2,i,j) = ((lambdal_relaxed + 2.d0/DIM*mul_relaxed)* &
        (de1(1,1,i,j) + de1(2,1,i,j)) + TWO * mul_relaxed * (de11(1,1,i,j) + de11(2,1,i,j))+ &
        (lambdalplus2mul_unrelaxed * (duxdx) + lambdal_unrelaxed* (duydy) ))

 dsigmayy(2,i,j) = ((lambdal_relaxed + 2.d0*mul_relaxed)* &
        (de1(1,1,i,j) + de1(2,1,i,j)) - TWO/DIM * mul_relaxed * (de11(1,1,i,j) + de11(2,1,i,j)) + &
        (lambdal_unrelaxed * (duxdx) + lambdalplus2mul_unrelaxed* (duydy) ))

! compute the stress using the unrelaxed Lame parameters (Carcione page 111)
  dsigmaxx_R(2,i,j) = lambdalplus2mul_relaxed * (duxdx) + lambdal_relaxed* (duydy)

  dsigmayy_R(2,i,j) = lambdal_relaxed * (duxdx) + lambdalplus2mul_relaxed* (duydy)

     enddo
    enddo

   do j=1,NY-1
     do i=2,NX
      mul_relaxed = 0.5d0 * (mu(i,j+1) + mu(i,j))
      mul_unrelaxed = mul_relaxed * Mu_nu2

        value_dvy_dx = ( c1 * (dvy(1,i,j) - dvy(1,i-1,j)) + c2 * (dvy(1,i+1,j) - dvy(1,i-2,j)) + &
            c3 * (dvy(1,i+2,j) - dvy(1,i-3,j)) + c4 * (dvy(1,i+3,j) - dvy(1,i-4,j)) )* ONE_OVER_DELTAX

         value_dvx_dy = ( c1 * (dvx(1,i,j+1) - dvx(1,i,j)) + c2 * (dvx(1,i,j+2) - dvx(1,i,j-1)) +  &
            c3 * (dvx(1,i,j+3) - dvx(1,i,j-2)) + c4 * (dvx(1,i,j+4) - dvx(1,i,j-3)) )* ONE_OVER_DELTAY

             duydx = value_dvy_dx
             duxdy = value_dvx_dy

           if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then
           memory_dvy_dx_1(2,i,j) = r_x_1(i) * memory_dvy_dx_1(1,i,j) + s_x_1(i) * value_dvy_dx
           memory_dvx_dy_1(2,i,j) = r_y_half_1(j) * memory_dvx_dy_1(1,i,j) + s_y_half_1(j) * value_dvx_dy

           duydx = value_dvy_dx / K_x_1(i)  + memory_dvy_dx_1(1,i,j)
           duxdy = value_dvx_dy / K_y_half_1(j) + memory_dvx_dy_1(1,i,j)
           endif

! evolution e12(1)
     tauinv = - inv_tau_sigma_nu2(1)
     Sn   = (duxdy+duydx) * phi_nu2(1)
     de12(1,2,i,j) = tauinv * de12(1,1,i,j) + Sn

! evolution e12(2)
     tauinv = - inv_tau_sigma_nu2(2)
     Sn   = (duxdy+duydx) * phi_nu2(2)
     de12(2,2,i,j) = tauinv * de12(2,1,i,j) + Sn

    dsigmaxy(2,i,j) = mul_relaxed * (de12(1,1,i,j) + de12(2,1,i,j))+mul_unrelaxed * (duxdy+duydx)
    dsigmaxy_R(2,i,j) = mul_relaxed * (duxdy+duydx)

      enddo
    enddo

        ! Update solution for next time step Delta_t
        vx(:,:) = vx(:,:) + rk42(inc) * dvx(2,:,:) * DELTAT
        vy(:,:) = vy(:,:) + rk42(inc) * dvy(2,:,:) * DELTAT
        sigmaxx(:,:) = sigmaxx(:,:) + rk42(inc) * dsigmaxx(2,:,:) * DELTAT
        sigmayy(:,:) = sigmayy(:,:) + rk42(inc) * dsigmayy(2,:,:) * DELTAT
        sigmaxy(:,:) = sigmaxy(:,:) + rk42(inc) * dsigmaxy(2,:,:) * DELTAT
        sigmaxx_R(:,:) = sigmaxx_R(:,:) + rk42(inc) * dsigmaxx_R(2,:,:) * DELTAT
        sigmayy_R(:,:) = sigmayy_R(:,:) + rk42(inc) * dsigmayy_R(2,:,:) * DELTAT
        sigmaxy_R(:,:) = sigmaxy_R(:,:) + rk42(inc) * dsigmaxy_R(2,:,:) * DELTAT

        e1(1,:,:) = e1(1,:,:) + rk42(inc) * de1(1,2,:,:) * DELTAT
        e1(2,:,:) = e1(2,:,:) + rk42(inc) * de1(2,2,:,:) * DELTAT
        e11(1,:,:) = e11(1,:,:) + rk42(inc) * de11(1,2,:,:) * DELTAT
        e11(2,:,:) = e11(2,:,:) + rk42(inc) * de11(2,2,:,:) * DELTAT
        e12(1,:,:) = e12(1,:,:) + rk42(inc) * de12(1,2,:,:) * DELTAT
        e12(2,:,:) = e12(2,:,:) + rk42(inc) * de12(2,2,:,:) * DELTAT

        memory_vx_dx_1(:,:) = memory_vx_dx_1(:,:) + rk42(inc) * memory_dvx_dx_1(2,:,:) * DELTAT
        memory_vx_dy_1(:,:) = memory_vx_dy_1(:,:) + rk42(inc) * memory_dvx_dy_1(2,:,:) * DELTAT
        memory_vy_dx_1(:,:) = memory_vy_dx_1(:,:) + rk42(inc) * memory_dvy_dx_1(2,:,:) * DELTAT
        memory_vy_dy_1(:,:) = memory_vy_dy_1(:,:) + rk42(inc) * memory_dvy_dy_1(2,:,:) * DELTAT
        memory_sigmaxx_dx_1(:,:) = memory_sigmaxx_dx_1(:,:) + rk42(inc) * memory_dsigmaxx_dx_1(2,:,:) * DELTAT
        memory_sigmayy_dy_1(:,:) = memory_sigmayy_dy_1(:,:) + rk42(inc) * memory_dsigmayy_dy_1(2,:,:) * DELTAT
        memory_sigmaxy_dx_1(:,:) = memory_sigmaxy_dx_1(:,:) + rk42(inc) * memory_dsigmaxy_dx_1(2,:,:) * DELTAT
        memory_sigmaxy_dy_1(:,:) = memory_sigmaxy_dy_1(:,:) + rk42(inc) * memory_dsigmaxy_dy_1(2,:,:) * DELTAT

        ! Dirichlet conditions (rigid boundaries) on all the edges of the grid
        dvx(:,-4:1,:) = ZERO
        dvx(:,NX:NX+4,:) = ZERO

        dvx(:,:,-4:1) = ZERO
        dvx(:,:,NY:NY+4) = ZERO

        dvy(:,-4:1,:) = ZERO
        dvy(:,NX:NX+4,:) = ZERO

        dvy(:,:,-4:1) = ZERO
        dvy(:,:,NY:NY+4) = ZERO

        vx(-4:1,:) = ZERO
        vx(:,-4:1) = ZERO
        vy(-4:1,:) = ZERO
        vy(:,-4:1) = ZERO

        vx(NX:NX+4,:) = ZERO
        vx(:,NY:NY+4) = ZERO
        vy(NX:NX+4,:) = ZERO
        vy(:,NY:NY+4) = ZERO

  enddo

  !vx(:,:) =  dvx(1,:,:)
  !vy(:,:) =  dvy(1,:,:)
  !sigmaxx(:,:) =  dsigmaxx(1,:,:)
  !sigmayy(:,:) =  dsigmayy(1,:,:)
  !sigmaxy(:,:) =  dsigmaxy(1,:,:)
  !sigmaxx_R(:,:) =  dsigmaxx_R(1,:,:)
  !sigmayy_R(:,:) =  dsigmayy_R(1,:,:)
  !sigmaxy_R(:,:) =  dsigmaxy_R(1,:,:)

  !e1(1,:,:) = de1(1,1,:,:)
  !e1(2,:,:) = de1(2,1,:,:)
  !e11(1,:,:) = de11(1,1,:,:)
  !e11(2,:,:) = de11(2,1,:,:)
  !e12(1,:,:) = de12(1,1,:,:)
  !e12(2,:,:) = de12(2,1,:,:)

! end of RK4 loop

! store seismograms
    do irec = 1,NREC
      sisvx(it,irec) = vx(ix_rec(irec),iy_rec(irec))
      sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec))
    enddo

! compute total energy in the medium (without the PML layers)
  local_energy_kinetic = ZERO
  local_energy_potential = ZERO


    do j = NPOINTS_PML, NY-NPOINTS_PML+1
      do i = NPOINTS_PML, NX-NPOINTS_PML+1

! compute kinetic energy first, defined as 1/2 rho ||v||^2
! in principle we should use rho_half_x_half_y instead of rho for vy
! in order to interpolate density at the right location in the staggered grid cell
! but in a homogeneous medium we can safely ignore it
      local_energy_kinetic = local_energy_kinetic + 0.5d0 * rho(i,j)*( &
              vx(i,j)**2 + vy(i,j)**2)

        total_energy_kinetic(it) = local_energy_kinetic

! add potential energy, defined as 1/2 epsilon_ij sigma_ij
! in principle we should interpolate the medium parameters at the right location
! in the staggered grid cell but in a homogeneous medium we can safely ignore it

! compute total field from split components
      epsilon_xx = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmaxx_R(i,j) - lambda(i,j) * sigmayy_R(i,j)) / &
                   (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)))

      epsilon_yy = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmayy_R(i,j) - lambda(i,j) * sigmaxx_R(i,j)) / &
                   (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)))

      epsilon_xy = sigmaxy_R(i,j) / (2.d0 * mu(i,j))

      local_energy_potential = local_energy_potential + &
        0.5d0 * (epsilon_xx * sigmaxx_R(i,j) + epsilon_yy * sigmayy_R(i,j) + &
        epsilon_yy * sigmayy_R(i,j)+ 2.d0 * epsilon_xy * sigmaxy_R(i,j))

      total_energy_potential(it) = local_energy_potential

        enddo
    enddo

      total_energy(it) = total_energy_kinetic(it) + total_energy_potential(it)

! output information
  if (mod(it,IT_DISPLAY) == 0 .or. it == 5) then
        Vsolidnorm = maxval(sqrt(vx**2 + vy**2))
      print *,'Time step # ',it,' out of ',NSTEP
      print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
      print *,'Max norm velocity vector V (m/s) = ',Vsolidnorm
      print *,'Total energy = ',total_energy(it)
! check stability of the code, exit if unstable
      if (Vsolidnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up in solid'

! save energy
     open(unit=21,file='energy.dat',status='unknown')
     do it2=1,NSTEP
       write(21,*) sngl(dble(it2-1)*DELTAT),total_energy_kinetic(it2), &
          total_energy_potential(it2),total_energy(it2)
     enddo
     close(21)

! save seismograms
    print *,'saving seismograms'
    print *
    call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT,t0)

    call create_color_image(vx(1:NX,1:NY),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1,max_amplitudeVx,JINTERFACE)
    call create_color_image(vy(1:NX,1:NY),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2,max_amplitudeVy,JINTERFACE)

    endif

! --- end of time loop
  enddo

! save seismograms
  call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT,t0)

! save total energy
  open(unit=20,file='RK4_energy.dat',status='unknown')
  do it = 1,NSTEP
    write(20,*) sngl(dble(it-1)*DELTAT),sngl(total_energy_kinetic(it)), &
            sngl(total_energy_potential(it)),sngl(total_energy(it))
  enddo
  close(20)

! create script for Gnuplot for total energy
  open(unit=20,file='RK4_plot_energy',status='unknown')
  write(20,*) 'set term x11'
  write(20,*) 'set term postscript landscape monochrome dashed "Helvetica" 22'
  write(20,*)
  write(20,*) 'set xlabel "Time (s)"'
  write(20,*) 'set ylabel "Total energy"'
  write(20,*)
  write(20,*) 'set output "ADEPML2D_total_energy_semilog.eps"'
  write(20,*) 'set logscale y'
  write(20,*) 'plot "RK4_energy.dat" t ''Total energy'' w l 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)
  close(20)

! create script for Gnuplot
  open(unit=20,file='plotgnu',status='unknown')
  write(20,*) 'set term x11'
  write(20,*) 'set term postscript landscape monochrome dashed "Helvetica" 22'
  write(20,*)
  write(20,*) 'set xlabel "Time (s)"'
  write(20,*) 'set ylabel "Amplitude (m / s)"'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vx_receiver_001.eps"'
  write(20,*) 'plot "RK4_Vx_file_001.dat" t ''Vx ADE-PML RK4'' w l 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_001.eps"'
  write(20,*) 'plot "RK4_Vy_file_001.dat" t ''Vy ADE-PML RK4'' w l 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vx_receiver_002.eps"'
  write(20,*) 'plot "RK4_Vx_file_002.dat" t ''Vx ADE-PML RK4'' w l 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_002.eps"'
  write(20,*) 'plot "RK4_Vy_file_002.dat" t ''Vy ADE-PML RK4'' w l 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vx_receiver_003.eps"'
  write(20,*) 'plot "RK4_Vx_file_003.dat" t ''Vx ADE-PML RK4'' w l 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_003.eps"'
  write(20,*) 'plot "RK4_Vy_file_003.dat" t ''Vy ADE-PML RK4'' w l 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  close(20)

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

!    write time stamp file to give information about progression of simulation
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

  print *
  print *,'End of the simulation'
  print *

  end program seismic_ADEPML_2D_viscoelastic_RK4_eighth_order

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
    write(file_name,"('RK4_Vx_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT-t0),' ',sngl(sisvx(it,irec))
    enddo
    close(11)
  enddo

! Y component
  do irec=1,nrec
    write(file_name,"('RK4_Vy_file_',i3.3,'.dat')") irec
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
              NPOINTS_PML,USE_PML_LEFT,USE_PML_RIGHT,USE_PML_BOTTOM,USE_PML_TOP,field_number,max_amplitude,JINTERFACE)


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
  logical USE_PML_LEFT,USE_PML_RIGHT,USE_PML_BOTTOM,USE_PML_TOP

  double precision, dimension(NX,NY) :: image_data_2D

  integer, dimension(nrec) :: ix_rec,iy_rec

  integer ix,iy,irec,JINTERFACE

  double precision max_amplitude

  character(len=100) file_name,system_command

  double precision normalized_value
  integer :: R, G, B

! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org
  if (field_number == 1) then
    write(file_name,"('image',i6.6,'_Vx.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vx.pnm image',i6.6,'_Vx.gif ; rm image',i6.6,'_Vx.pnm')") it,it,it
  endif
  if (field_number == 2) then
    write(file_name,"('image',i6.6,'_Vy.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vy.pnm image',i6.6,'_Vy.gif ; rm image',i6.6,'_Vy.pnm')") it,it,it
  endif
  if (field_number == 3) then
    write(file_name,"('image',i6.6,'_Vnorm.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vnorm.pnm image',i6.6,'_Vnorm.gif ; rm image',i6.6,'_Vnorm.pnm')") it,it,it
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
  else if ((USE_PML_LEFT .and. ix == NPOINTS_PML) .or. &
          (USE_PML_RIGHT .and. ix == NX - NPOINTS_PML) .or. &
          (USE_PML_BOTTOM .and. iy == NPOINTS_PML) .or. &
          (USE_PML_TOP .and. iy == NY - NPOINTS_PML)) then
      R = 255
      G = 150
      B = 0
 else if (iy == JINTERFACE) then
        R = 0
        G = 0
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

! call the system to convert image to JPEG
! call system(system_command)

  end subroutine create_color_image

