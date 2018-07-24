!
! SEISMIC_CPML Version 1.1.3, July 2018.
!
! Copyright CNRS, France.
! Contributor: Dimitri Komatitsch, komatitsch aT lma DOT cnrs-mrs DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional heterogeneous isotropic viscoelastic wave equation
! using a finite-difference method with Convolutional Perfectly Matched
! Layer (C-PML) conditions.
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

  program seismic_CPML_2D_viscoelast_fourth

! 2D finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions for an heterogeneous isotropic viscoelastic medium

! Dimitri Komatitsch, CNRS, Marseille, July 2018.

! A fourth-order spatially-staggered grid formulation is used:
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
!        e13 |         |         |
!    (memory |         |         |
!  variable) |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!           v_x    sigma_xx
!                  sigma_yy
!                  e1 (viscoelastic memory variable)
!                  e11 (viscoelastic memory variable)
!

! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
! If you use this code for your own research, please cite some (or all) of these
! articles:
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

! include viscoelastic attenuation or not
  logical, parameter :: VISCOELASTIC_ATTENUATION = .true.

! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 2001
  integer, parameter :: NY = 2001

! size of a grid cell
  double precision, parameter :: DELTAX = 1.5d0
  double precision, parameter :: DELTAY = DELTAX

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 10

! P-velocity and density
! the unrelaxed value is the value at frequency = 0 (the relaxed value would be the value at frequency = +infinity)
  double precision, parameter :: cp_unrelaxed = 2000.d0
  double precision, parameter :: cs_unrelaxed = cp_unrelaxed / 1.732d0
  double precision, parameter :: density = 2000.d0

! Time step in seconds.
! The CFL stability number for the O(2,2) algorithm is 1 / sqrt(2) = 0.707
! i.e. one must choose  cp * deltat / deltax < 0.707.
! For the O(2,4) algorithm used here it is a bit more restrictive,
! it is cp * deltat / deltax < 0.606  (see Levander 1988 eq (7)).
! However this only ensures that the scheme is stable. To have a scheme that is both stable and accurate,
! for O(2,4) some numerical tests show that one needs to take about half of that,
! i.e. choose deltat so that cp * deltat / deltax is equal to about 0.30 or so. (or any value below; but not above).
! Since the time scheme is only second order, this also depends on how many time steps are performed in total
! (i.e. what the value of NSTEP below is); for large values of NSTEP, of course numerical errors will start to accumulate.
  double precision, parameter :: DELTAT = 2.2d-4

! total number of time steps
  integer, parameter :: NSTEP = 5200

! parameters for the source
  double precision, parameter :: f0 = 35.d0
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.d0

! source (force)
  double precision, parameter :: xsource = 1500.d0
  double precision, parameter :: ysource = 1500.d0
  integer, parameter :: ISOURCE = xsource / DELTAX + 1
  integer, parameter :: JSOURCE = ysource / DELTAY + 1
! angle of source force in degrees and clockwise, with respect to the vertical (Y) axis
  double precision, parameter :: ANGLE_FORCE = 0.d0

! receivers
  integer, parameter :: NREC = 1
!! DK DK I use 2301 here instead of 2300 in order to fall exactly on a grid point
  double precision, parameter :: xdeb = 2301.d0   ! first receiver x in meters
  double precision, parameter :: ydeb = 2301.d0   ! first receiver y in meters
  double precision, parameter :: xfin = 2301.d0   ! last receiver x in meters
  double precision, parameter :: yfin = 2301.d0   ! last receiver y in meters

! to compute energy curves for the whole medium (optional, but useful e.g. to produce
! energy variation figures for articles); but expensive option, thus off by default
  logical, parameter :: COMPUTE_ENERGY = .false.

! display information on the screen from time to time
  integer, parameter :: IT_DISPLAY = 200

! compute some constants once and for all for the fourth-order spatial scheme
! These coefficients are given for instance by Levander, Geophysics, vol. 53(11), p. 1436, equation (A-2)
  double precision, parameter :: NINE_OVER_8_DELTAX = 9.d0 / (8.d0*DELTAX)
  double precision, parameter :: NINE_OVER_8_DELTAY = 9.d0 / (8.d0*DELTAY)
  double precision, parameter :: ONE_OVER_24_DELTAX = 1.d0 / (24.d0*DELTAX)
  double precision, parameter :: ONE_OVER_24_DELTAY = 1.d0 / (24.d0*DELTAY)

! value of PI
  double precision, parameter :: PI = 3.141592653589793238462643d0

! conversion from degrees to radians
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! zero
  double precision, parameter :: ZERO = 0.d0

  double precision, parameter :: TWO_THIRDS = 2.d0 / 3.d0

! large value for maximum
  double precision, parameter :: HUGEVAL = 1.d+30

! threshold above which we consider that the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
! in order to be able to use a fourth-order spatial operator on the edges of the model
! here we define the arrays with size (0:NX+1,0:NY+1) instead of size (NX,NY) as in the second-order case
  double precision, dimension(0:NX+1,0:NY+1) :: vx,vy,sigma_xx,sigma_yy,sigma_xy,lambda_unrelaxed,mu_unrelaxed,rho

! to interpolate material parameters or velocity at the right location in the staggered grid cell
  double precision :: lambda_half_x,mu_half_x,lambda_plus_mu_half_x,lambda_plus_two_mu_half_x,mu_half_y
  double precision :: rho_half_x_half_y,vy_interpolated

! for evolution of total energy in the medium
  double precision :: epsilon_xx,epsilon_yy,epsilon_xy
  double precision, dimension(NSTEP) :: total_energy_kinetic,total_energy_potential

! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
  double precision, parameter :: K_MAX_PML = 1.d0
  double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
  double precision, dimension(NX,NY) :: &
      memory_dvx_dx, &
      memory_dvx_dy, &
      memory_dvy_dx, &
      memory_dvy_dy, &
      memory_dsigma_xx_dx, &
      memory_dsigma_yy_dy, &
      memory_dsigma_xy_dx, &
      memory_dsigma_xy_dy

  double precision :: &
      value_dvx_dx, &
      value_dvx_dy, &
      value_dvy_dx, &
      value_dvy_dy, &
      value_dsigma_xx_dx, &
      value_dsigma_yy_dy, &
      value_dsigma_xy_dx, &
      value_dsigma_xy_dy

! 1D arrays for the damping profiles
  double precision, dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half, &
                                     one_over_K_x,one_over_K_x_half
  double precision, dimension(NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half, &
                                     one_over_K_y,one_over_K_y_half

  double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
  double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

! for the source
  double precision :: a,t,force_x,force_y,force_source_term

! for receivers
! Please note something important: the two components of the velocity vector are not defined at the same location,
! Vy is half a grid cell away from Vx (see ASCII figure at the beginning of this program).
! Thus this means there are "two receivers" rather than one, one recording Vx and another one, half a grid cell away, recording Vy.
! If you need to use both components in real applications (and of course we will),
! you will need to interpolate Vy to the location of Vx using:
!
! interpolate vy back at the location of vx, to be able to use both at the same location
!       vy_interpolated = 0.25d0 * (vy(i,j) + vy(i-1,j) + vy(i-1,j-1) + vy(i,j-1))
!
  double precision xspacerec,yspacerec,distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec
  double precision, dimension(NREC) :: xrec,yrec
  integer :: myNREC

! for seismograms
  double precision, dimension(NSTEP,NREC) :: sisvx,sisvy,sispressure

  integer :: i,j,it,irec

  double precision :: Courant_number,velocnorm

! for attenuation (viscoelasticity)

! attenuation quality factors Qp and Qs to use
! BEWARE: we use Qp and Qs here, not QKappa and Qmu.
! BEWARE: While Qmu is always equal to Qs, QKappa is not equal to Qp,
! BEWARE: to convert from one to the other if your input data have Qkappa and Qmu you can use
! BEWARE: the program conversion_between_Qp_Qs_and_Qkappa_Qmu_from_Dahlen_Tromp_959_960_in_3D_and_in_2D_plane_strain.f90
! BEWARE: that is included in this software package.
  double precision, parameter :: Qp = 65.d0
  double precision, parameter :: Qs = 55.d0

! number of Zener standard linear solids in parallel
  integer, parameter :: N_SLS = 3

! attenuation constants
  double precision, dimension(N_SLS) :: tau_epsilon_nu1,tau_sigma_nu1,one_over_tau_sigma_nu1, &
                           HALF_DELTAT_over_tau_sigma_nu1,multiplication_factor_tau_sigma_nu1,DELTAT_phi_nu1
  double precision, dimension(N_SLS) :: tau_epsilon_nu2,tau_sigma_nu2,one_over_tau_sigma_nu2, &
                           HALF_DELTAT_over_tau_sigma_nu2,multiplication_factor_tau_sigma_nu2,DELTAT_phi_nu2

! memory variable and other arrays for attenuation
  double precision, dimension(NX,NY,N_SLS) :: memory_variable_R_e1_dot,memory_variable_R_e1_dot_old
  double precision, dimension(NX,NY,N_SLS) :: memory_variable_R_e11_dot,memory_variable_R_e11_dot_old
  double precision, dimension(NX,NY,N_SLS) :: memory_variable_R_e13_dot,memory_variable_R_e13_dot_old
  integer :: i_sls
  double precision :: sum_of_memory_variables_e1,sum_of_memory_variables_e11,sum_of_memory_variables_e13

! this defines the typical frequency range in which we use optimization to find the tau values that fit a given Q in that band
  double precision :: f_min_attenuation,f_max_attenuation

!---
!--- program starts here
!---

  print *
  print *,'2D viscoelastic finite-difference code in velocity and stress formulation with C-PML'
  print *

! display size of the model
  print *
  print *,'NX = ',NX
  print *,'NY = ',NY
  print *
  print *,'size of the model along X = ',(NX - 1) * DELTAX
  print *,'size of the model along Y = ',(NY - 1) * DELTAY
  print *
  print *,'Total number of grid points = ',NX * NY
  print *

! for attenuation (viscoelasticity)
  if (VISCOELASTIC_ATTENUATION) then

  print *,'Qp quality factor used for attenuation = ',Qp
  print *,'Qs quality factor used for attenuation = ',Qs
  print *,'Number of Zener standard linear solids used to mimic the viscoelastic behavior (N_SLS) = ',N_SLS
  print *

! this defines the typical frequency range in which we use optimization to find the tau values that fit a given Q in that band
! f_min and f_max are computed as : f_max/f_min=12 and (log(f_min)+log(f_max))/2 = log(f0)
  f_min_attenuation = exp(log(f0)-log(12.d0)/2.d0)
  f_max_attenuation = 12.d0 * f_min_attenuation

! call the SolvOpt() nonlinear optimization routine to compute the tau_epsilon and tau_sigma values from a given Q factor
  print *,'Values for Qp:'
  print *
  call compute_attenuation_coeffs(N_SLS,Qp,f0,f_min_attenuation,f_max_attenuation,tau_epsilon_nu1,tau_sigma_nu1)
  print *,'Values for Qs:'
  print *
  call compute_attenuation_coeffs(N_SLS,Qs,f0,f_min_attenuation,f_max_attenuation,tau_epsilon_nu2,tau_sigma_nu2)

  else

! dummy values in the non-dissipative case
    tau_epsilon_nu1(:) = 1.d0
    tau_sigma_nu1(:) = 1.d0

    tau_epsilon_nu2(:) = 1.d0
    tau_sigma_nu2(:) = 1.d0

  endif

! precompute the inverse once and for all, to save computation time in the time loop below
! (on computers, a multiplication is very significantly cheaper than a division)
  one_over_tau_sigma_nu1(:) = 1.d0 / tau_sigma_nu1(:)
  one_over_tau_sigma_nu2(:) = 1.d0 / tau_sigma_nu2(:)

  HALF_DELTAT_over_tau_sigma_nu1(:) = 0.5d0 * DELTAT / tau_sigma_nu1(:)
  HALF_DELTAT_over_tau_sigma_nu2(:) = 0.5d0 * DELTAT / tau_sigma_nu2(:)

  multiplication_factor_tau_sigma_nu1(:) = 1.d0 / (1.d0 + 0.5d0 * DELTAT * one_over_tau_sigma_nu1(:))
  multiplication_factor_tau_sigma_nu2(:) = 1.d0 / (1.d0 + 0.5d0 * DELTAT * one_over_tau_sigma_nu2(:))

  ! use the right formula with 1/N included
  DELTAT_phi_nu1(:) = DELTAT * (1.d0 - tau_epsilon_nu1(:)/tau_sigma_nu1(:)) / tau_sigma_nu1(:) / sum(tau_epsilon_nu1/tau_sigma_nu1)
  DELTAT_phi_nu2(:) = DELTAT * (1.d0 - tau_epsilon_nu2(:)/tau_sigma_nu2(:)) / tau_sigma_nu2(:) / sum(tau_epsilon_nu2/tau_sigma_nu2)

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
  thickness_PML_x = NPOINTS_PML * DELTAX
  thickness_PML_y = NPOINTS_PML * DELTAY

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0

! check that NPOWER is okay
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * cp_unrelaxed * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * cp_unrelaxed * log(Rcoef) / (2.d0 * thickness_PML_y)

  print *,'d0_x = ',d0_x
  print *,'d0_y = ',d0_y
  print *

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

! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x
  xoriginright = (NX-1)*DELTAX - thickness_PML_x

  do i = 1,NX

! abscissa of current grid point along the damping profile
    xval = DELTAX * dble(i-1)

!---------- left edge
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

!---------- right edge
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

!---------- bottom edge
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

!---------- top edge
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

! precompute the inverse once and for all, to save computation time in the time loop below
! (on computers, a multiplication is very significantly cheaper than a division)
  one_over_K_x(:) = 1.d0 / K_x(:)
  one_over_K_x_half(:) = 1.d0 / K_x_half(:)
  one_over_K_y(:) = 1.d0 / K_y(:)
  one_over_K_y_half(:) = 1.d0 / K_y_half(:)

! compute the Lame parameter and density
  do j = 1,NY
    do i = 1,NX
      rho(i,j) = density
      mu_unrelaxed(i,j) = density*cs_unrelaxed*cs_unrelaxed
      lambda_unrelaxed(i,j) = density*cp_unrelaxed*cp_unrelaxed - 2.d0*mu_unrelaxed(i,j)
    enddo
  enddo

! print position of the source
  print *,'Position of the source:'
  print *
  print *,'x = ',xsource
  print *,'y = ',ysource
  print *

! define location of receivers
  print *,'There are ',nrec,' receivers'
  print *
  if (NREC > 1) then
! this is to avoid a warning with GNU gfortran at compile time about division by zero when NREC = 1
    myNREC = NREC
    xspacerec = (xfin-xdeb) / dble(myNREC-1)
    yspacerec = (yfin-ydeb) / dble(myNREC-1)
  else
    xspacerec = 0.d0
    yspacerec = 0.d0
  endif
  do irec=1,nrec
    xrec(irec) = xdeb + dble(irec-1)*xspacerec
    yrec(irec) = ydeb + dble(irec-1)*yspacerec
  enddo

! find closest grid point for each receiver
  do irec=1,nrec
    dist = HUGEVAL
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((DELTAX*dble(i-1) - xrec(irec))**2 + (DELTAY*dble(j-1) - yrec(irec))**2)
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
! R. Courant, K. O. Friedrichs and H. Lewy (1928)
! For this O(2,4) scheme, when DELTAX == DELTAY the Courant number is given by Levander, Geophysics, vol. 53(11), p. 1427,
! equation (7) and is equal to 0.606 (it is thus smaller than that of the O(2,2) scheme, which is 1/sqrt(2) = 0.707,
! i.e. when switching to a fourth-order spatial scheme one needs a time step that is about 0.707 / 0.606 = 1.167 times smaller.
  if (DELTAX == DELTAY) then
    Courant_number = cp_unrelaxed * DELTAT / DELTAX
    print *,'Courant number is ',Courant_number
    print *,' (the maximum possible value is 0.606; in practice for accuracy reasons a value not larger than 0.30 is recommended)'
    print *
    if (Courant_number > 0.606) stop 'time step is too large, simulation will be unstable'
  endif

! suppress old files (can be commented out if "call system" is missing in your compiler)
  call system('rm -f Vx_file*.dat Vy_file*.dat image*.pnm image*.gif')

! initialize arrays
  vx(:,:) = ZERO
  vy(:,:) = ZERO
  sigma_xx(:,:) = ZERO
  sigma_yy(:,:) = ZERO
  sigma_xy(:,:) = ZERO
  memory_variable_R_e1_dot(:,:,:) = ZERO
  memory_variable_R_e1_dot_old(:,:,:) = ZERO
  memory_variable_R_e11_dot(:,:,:) = ZERO
  memory_variable_R_e11_dot_old(:,:,:) = ZERO
  memory_variable_R_e13_dot(:,:,:) = ZERO
  memory_variable_R_e13_dot_old(:,:,:) = ZERO

! PML
  memory_dvx_dx(:,:) = ZERO
  memory_dvx_dy(:,:) = ZERO
  memory_dvy_dx(:,:) = ZERO
  memory_dvy_dy(:,:) = ZERO
  memory_dsigma_xx_dx(:,:) = ZERO
  memory_dsigma_yy_dy(:,:) = ZERO
  memory_dsigma_xy_dx(:,:) = ZERO
  memory_dsigma_xy_dy(:,:) = ZERO

! initialize seismograms
  sisvx(:,:) = ZERO
  sisvy(:,:) = ZERO
  sispressure(:,:) = ZERO

! initialize total energy
  total_energy_kinetic(:) = ZERO
  total_energy_potential(:) = ZERO

  if (VISCOELASTIC_ATTENUATION) then
    print *,'adding VISCOELASTIC_ATTENUATION (i.e., running a viscoelastic simulation)'
  else
    print *,'not adding VISCOELASTIC_ATTENUATION (i.e., running a purely elastic simulation)'
  endif
  print *

!---
!---  beginning of time loop
!---

  do it = 1,NSTEP

!-----------------------------------------------------------------------
! compute the stress tensor and update memory variables for C-PML
! also update memory variables for viscoelastic attenuation if needed
!-----------------------------------------------------------------------

! we purposely leave this "if" test outside of the loops to make sure the compiler can optimize these loops;
! with an "if" test inside most compilers cannot
  if (.not. VISCOELASTIC_ATTENUATION) then

    do j = 2,NY
      do i = 1,NX-1

! interpolate material parameters at the right location in the staggered grid cell
      lambda_half_x = 0.5d0 * (lambda_unrelaxed(i+1,j) + lambda_unrelaxed(i,j))
      mu_half_x = 0.5d0 * (mu_unrelaxed(i+1,j) + mu_unrelaxed(i,j))
      lambda_plus_two_mu_half_x = lambda_half_x + 2.d0 * mu_half_x

      value_dvx_dx = (vx(i+1,j) - vx(i,j)) * NINE_OVER_8_DELTAX + (vx(i-1,j) - vx(i+2,j)) * ONE_OVER_24_DELTAX
      value_dvy_dy = (vy(i,j) - vy(i,j-1)) * NINE_OVER_8_DELTAY + (vy(i,j-2) - vy(i,j+1)) * ONE_OVER_24_DELTAY

      memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx
      memory_dvy_dy(i,j) = b_y(j) * memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy

      value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
      value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j)

      sigma_xx(i,j) = sigma_xx(i,j) + (lambda_plus_two_mu_half_x * value_dvx_dx + lambda_half_x * value_dvy_dy) * DELTAT

      sigma_yy(i,j) = sigma_yy(i,j) + (lambda_half_x * value_dvx_dx + lambda_plus_two_mu_half_x * value_dvy_dy) * DELTAT

      enddo
    enddo

    do j = 1,NY-1
      do i = 2,NX

! interpolate material parameters at the right location in the staggered grid cell
        mu_half_y = 0.5d0 * (mu_unrelaxed(i,j+1) + mu_unrelaxed(i,j))

        value_dvy_dx = (vy(i,j) - vy(i-1,j)) * NINE_OVER_8_DELTAX + (vy(i-2,j) - vy(i+1,j)) * ONE_OVER_24_DELTAX
        value_dvx_dy = (vx(i,j+1) - vx(i,j)) * NINE_OVER_8_DELTAY + (vx(i,j-1) - vx(i,j+2)) * ONE_OVER_24_DELTAY

        memory_dvy_dx(i,j) = b_x(i) * memory_dvy_dx(i,j) + a_x(i) * value_dvy_dx
        memory_dvx_dy(i,j) = b_y_half(j) * memory_dvx_dy(i,j) + a_y_half(j) * value_dvx_dy

        value_dvy_dx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j)
        value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j)

        sigma_xy(i,j) = sigma_xy(i,j) + mu_half_y * (value_dvy_dx + value_dvx_dy) * DELTAT

      enddo
    enddo

  else

! the present becomes the past for the memory variables.
! in C or C++ we could replace this with an exchange of pointers on the arrays
! in order to avoid a memory copy of the whole array.
    memory_variable_R_e1_dot_old(:,:,:) = memory_variable_R_e1_dot(:,:,:)
    memory_variable_R_e11_dot_old(:,:,:) = memory_variable_R_e11_dot(:,:,:)
    memory_variable_R_e13_dot_old(:,:,:) = memory_variable_R_e13_dot(:,:,:)

    do j = 2,NY
      do i = 1,NX-1

! interpolate material parameters at the right location in the staggered grid cell
        lambda_half_x = 0.5d0 * (lambda_unrelaxed(i+1,j) + lambda_unrelaxed(i,j))
        mu_half_x = 0.5d0 * (mu_unrelaxed(i+1,j) + mu_unrelaxed(i,j))
        lambda_plus_mu_half_x = lambda_half_x + mu_half_x
        lambda_plus_two_mu_half_x = lambda_half_x + 2.d0 * mu_half_x

        value_dvx_dx = (vx(i+1,j) - vx(i,j)) * NINE_OVER_8_DELTAX + (vx(i-1,j) - vx(i+2,j)) * ONE_OVER_24_DELTAX
        value_dvy_dy = (vy(i,j) - vy(i,j-1)) * NINE_OVER_8_DELTAY + (vy(i,j-2) - vy(i,j+1)) * ONE_OVER_24_DELTAY

        memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx
        memory_dvy_dy(i,j) = b_y(j) * memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy

        value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
        value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j)

! use the Auxiliary Differential Equation form, which is second-order accurate in time if implemented following
! eq (14) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994), which is what we do here
        sum_of_memory_variables_e1 = 0.d0
        sum_of_memory_variables_e11 = 0.d0
        do i_sls = 1,N_SLS
! this average of the two terms comes from eq (14) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
          memory_variable_R_e1_dot(i,j,i_sls) = (memory_variable_R_e1_dot_old(i,j,i_sls) + &
                   (value_dvx_dx + value_dvy_dy) * DELTAT_phi_nu1(i_sls) - &
                   memory_variable_R_e1_dot_old(i,j,i_sls) * HALF_DELTAT_over_tau_sigma_nu1(i_sls)) &
                      * multiplication_factor_tau_sigma_nu1(i_sls)

          memory_variable_R_e11_dot(i,j,i_sls) = (memory_variable_R_e11_dot_old(i,j,i_sls) + &
                   0.5d0 * (value_dvx_dx - value_dvy_dy) * DELTAT_phi_nu2(i_sls) - &
                   memory_variable_R_e11_dot_old(i,j,i_sls) * HALF_DELTAT_over_tau_sigma_nu2(i_sls)) &
                      * multiplication_factor_tau_sigma_nu2(i_sls)

          sum_of_memory_variables_e1 = sum_of_memory_variables_e1 + &
                      memory_variable_R_e1_dot(i,j,i_sls) + memory_variable_R_e1_dot_old(i,j,i_sls)

          sum_of_memory_variables_e11 = sum_of_memory_variables_e11 + &
                      memory_variable_R_e11_dot(i,j,i_sls) + memory_variable_R_e11_dot_old(i,j,i_sls)
        enddo

        sigma_xx(i,j) = sigma_xx(i,j) + &
           (lambda_plus_two_mu_half_x * value_dvx_dx + lambda_half_x * value_dvy_dy &
! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
! this average of the two terms comes from eq (13) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
          + (0.5d0 * lambda_plus_mu_half_x * sum_of_memory_variables_e1 + mu_half_x * sum_of_memory_variables_e11)) * DELTAT

        sigma_yy(i,j) = sigma_yy(i,j) + &
           (lambda_half_x * value_dvx_dx + lambda_plus_two_mu_half_x * value_dvy_dy &
! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
! this average of the two terms comes from eq (13) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
          + (0.5d0 * lambda_plus_mu_half_x * sum_of_memory_variables_e1 - mu_half_x * sum_of_memory_variables_e11)) * DELTAT

      enddo
    enddo

    do j = 1,NY-1
      do i = 2,NX

! interpolate material parameters at the right location in the staggered grid cell
        mu_half_y = 0.5d0 * (mu_unrelaxed(i,j+1) + mu_unrelaxed(i,j))

        value_dvy_dx = (vy(i,j) - vy(i-1,j)) * NINE_OVER_8_DELTAX + (vy(i-2,j) - vy(i+1,j)) * ONE_OVER_24_DELTAX
        value_dvx_dy = (vx(i,j+1) - vx(i,j)) * NINE_OVER_8_DELTAY + (vx(i,j-1) - vx(i,j+2)) * ONE_OVER_24_DELTAY

        memory_dvy_dx(i,j) = b_x(i) * memory_dvy_dx(i,j) + a_x(i) * value_dvy_dx
        memory_dvx_dy(i,j) = b_y_half(j) * memory_dvx_dy(i,j) + a_y_half(j) * value_dvx_dy

        value_dvy_dx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j)
        value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j)

! use the Auxiliary Differential Equation form, which is second-order accurate in time if implemented following
! eq (14) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994), which is what we do here
        sum_of_memory_variables_e13 = 0.d0
        do i_sls = 1,N_SLS
! this average of the two terms comes from eq (14) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
          memory_variable_R_e13_dot(i,j,i_sls) = (memory_variable_R_e13_dot_old(i,j,i_sls) + &
                   (value_dvy_dx + value_dvx_dy) * DELTAT_phi_nu2(i_sls) - &
                   memory_variable_R_e13_dot_old(i,j,i_sls) * HALF_DELTAT_over_tau_sigma_nu2(i_sls)) &
                      * multiplication_factor_tau_sigma_nu2(i_sls)

          sum_of_memory_variables_e13 = sum_of_memory_variables_e13 + &
                      memory_variable_R_e13_dot(i,j,i_sls) + memory_variable_R_e13_dot_old(i,j,i_sls)
        enddo

        sigma_xy(i,j) = sigma_xy(i,j) + mu_half_y * (value_dvy_dx + value_dvx_dy &
! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
! this average of the two terms comes from eq (13) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
                                      + 0.5d0 * sum_of_memory_variables_e13) * DELTAT

      enddo
    enddo

  endif

!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

  do j = 2,NY
    do i = 2,NX

      value_dsigma_xx_dx = (sigma_xx(i,j) - sigma_xx(i-1,j)) * NINE_OVER_8_DELTAX + &
                                   (sigma_xx(i-2,j) - sigma_xx(i+1,j)) * ONE_OVER_24_DELTAX
      value_dsigma_xy_dy = (sigma_xy(i,j) - sigma_xy(i,j-1)) * NINE_OVER_8_DELTAY + &
                                   (sigma_xy(i,j-2) - sigma_xy(i,j+1)) * ONE_OVER_24_DELTAY

      memory_dsigma_xx_dx(i,j) = b_x(i) * memory_dsigma_xx_dx(i,j) + a_x(i) * value_dsigma_xx_dx
      memory_dsigma_xy_dy(i,j) = b_y(j) * memory_dsigma_xy_dy(i,j) + a_y(j) * value_dsigma_xy_dy

      value_dsigma_xx_dx = value_dsigma_xx_dx / K_x(i) + memory_dsigma_xx_dx(i,j)
      value_dsigma_xy_dy = value_dsigma_xy_dy / K_y(j) + memory_dsigma_xy_dy(i,j)

      vx(i,j) = vx(i,j) + (value_dsigma_xx_dx + value_dsigma_xy_dy) * DELTAT / rho(i,j)

    enddo
  enddo

  do j = 1,NY-1
    do i = 1,NX-1

! interpolate density at the right location in the staggered grid cell
      rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))

      value_dsigma_xy_dx = (sigma_xy(i+1,j) - sigma_xy(i,j)) * NINE_OVER_8_DELTAX + &
                                   (sigma_xy(i-1,j) - sigma_xy(i+2,j)) * ONE_OVER_24_DELTAX

      value_dsigma_yy_dy = (sigma_yy(i,j+1) - sigma_yy(i,j)) * NINE_OVER_8_DELTAY + &
                                   (sigma_yy(i,j-1) - sigma_yy(i,j+2)) * ONE_OVER_24_DELTAY

      memory_dsigma_xy_dx(i,j) = b_x_half(i) * memory_dsigma_xy_dx(i,j) + a_x_half(i) * value_dsigma_xy_dx
      memory_dsigma_yy_dy(i,j) = b_y_half(j) * memory_dsigma_yy_dy(i,j) + a_y_half(j) * value_dsigma_yy_dy

      value_dsigma_xy_dx = value_dsigma_xy_dx / K_x_half(i) + memory_dsigma_xy_dx(i,j)
      value_dsigma_yy_dy = value_dsigma_yy_dy / K_y_half(j) + memory_dsigma_yy_dy(i,j)

      vy(i,j) = vy(i,j) + (value_dsigma_xy_dx + value_dsigma_yy_dy) * DELTAT / rho_half_x_half_y

    enddo
  enddo

! add the source (force vector located at a given grid point)
  a = pi*pi*f0*f0
  t = dble(it-1)*DELTAT

! Gaussian
! force_source_term = - factor * exp(-a*(t-t0)**2) / (2.d0 * a)

! first derivative of a Gaussian
! force_source_term = factor * (t-t0)*exp(-a*(t-t0)**2)

! Ricker source time function (second derivative of a Gaussian)
  force_source_term = factor * (1.d0 - 2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)

! to get the right amplitude of the force, we need to divide by the area of a grid cell
! (we checked that against the analytical solution in a homogeneous medium for a force source)
  force_source_term = force_source_term / (DELTAX * DELTAY)

! define location of the source
  i = ISOURCE
  j = JSOURCE

  force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * force_source_term
  force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * force_source_term

! interpolate density at the right location in the staggered grid cell
  rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))

! we want seismograms to be representing velocity, for the case of the seismic wave equation
! representing displacement for a Ricker (i.e., second derivative of a Gaussian) source in displacement.
! Since the force source is added to d(velocity)/dt in this split velocity and stress scheme
! we need to select the second derivative of a Gaussian as a source time wavelet
! by analogy with a Ricker (i.e. a second derivative) added to d2(displacement)/dt2
! as in the unsplit equation written in displacement only.
! Since the formula is d(velocity)/dt = (velocity_new - velocity_old) / DELTAT = force_source_term
! we also need to multiply by DELTAT here to avoid having an amplitude of the seismogram
! that varies when one changes the time step, i.e. we write:
! velocity_new = velocity_old + force_source_term * DELTAT at the source grid point
  vx(i,j) = vx(i,j) + force_x * DELTAT / rho(i,j)
  vy(i,j) = vy(i,j) + force_y * DELTAT / rho_half_x_half_y

! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  vx(1,:) = ZERO
  vx(NX,:) = ZERO

  vx(:,1) = ZERO
  vx(:,NY) = ZERO

  vy(1,:) = ZERO
  vy(NX,:) = ZERO

  vy(:,1) = ZERO
  vy(:,NY) = ZERO

! store seismograms
  do irec = 1,NREC

! beware here that the two components of the velocity vector are not defined at the same point
! in a staggered grid, and thus the two components of the velocity vector are recorded at slightly different locations,
! vy is staggered by half a grid cell along X and along Y with respect to vx
    sisvx(it,irec) = vx(ix_rec(irec),iy_rec(irec))
    sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec))

! from L. S. Bennethum, Compressibility Moduli for Porous Materials Incorporating Volume Fraction,
! J. Engrg. Mech., vol. 132(11), p. 1205-1214 (2006), below equation (5):
! for a 3D isotropic solid, pressure is defined in terms of the trace of the stress tensor as
! p = -1/3 (t11 + t22 + t33) where t is the Cauchy stress tensor.

! to compute pressure in 3D in an elastic solid, one uses pressure = - trace(sigma) / 3
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_yy
! sigma_zz = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_zz
! pressure = - trace(sigma) / 3 = - (lambda + 2/3 mu) trace(epsilon) = - kappa * trace(epsilon)
!
! to compute pressure in 2D in an elastic solid in the plane strain convention i.e. in the P-SV case,
! one still uses pressure = - trace(sigma) / 3 but taking into account the fact
! that the off-plane strain epsilon_zz is zero by definition of the plane strain convention
! but thus the off-plane stress sigma_zz is not equal to zero,
! one has instead:  sigma_zz = lambda * (epsilon_xx + epsilon_yy), thus
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_yy
! sigma_zz = lambda * (epsilon_xx + epsilon_yy)
! pressure = - trace(sigma) / 3 = - (lambda + 2*mu/3) (epsilon_xx + epsilon_yy)

    i = ix_rec(irec)
    j = iy_rec(irec)

! interpolate material parameters at the right location in the staggered grid cell
    lambda_half_x = 0.5d0 * (lambda_unrelaxed(i+1,j) + lambda_unrelaxed(i,j))
    mu_half_x = 0.5d0 * (mu_unrelaxed(i+1,j) + mu_unrelaxed(i,j))
    epsilon_xx = ((lambda_half_x + 2.d0*mu_half_x) * sigma_xx(i,j) - lambda_half_x * &
      sigma_yy(i,j)) / (4.d0 * mu_half_x * (lambda_half_x + mu_half_x))
    epsilon_yy = ((lambda_half_x + 2.d0*mu_half_x) * sigma_yy(i,j) - lambda_half_x * &
      sigma_xx(i,j)) / (4.d0 * mu_half_x * (lambda_half_x + mu_half_x))

    sispressure(it,irec) = - (lambda_half_x + TWO_THIRDS*mu_half_x) * (epsilon_xx + epsilon_yy)

  enddo

! compute total energy in the medium (without the PML layers)
  if (COMPUTE_ENERGY) then

! compute kinetic energy first, defined as 1/2 rho ||v||^2
    total_energy_kinetic(it) = ZERO
    do j = NPOINTS_PML+1, NY-NPOINTS_PML
      do i = NPOINTS_PML+1, NX-NPOINTS_PML
! interpolate vy back at the location of vx, to be able to use both at the same location
        vy_interpolated = 0.25d0 * (vy(i,j) + vy(i-1,j) + vy(i-1,j-1) + vy(i,j-1))
        total_energy_kinetic(it) = total_energy_kinetic(it) + 0.5d0 * rho(i,j) * (vx(i,j)**2 + vy_interpolated**2)
      enddo
    enddo

! add potential energy, defined as 1/2 epsilon_ij sigma_ij
    total_energy_potential(it) = ZERO
    do j = NPOINTS_PML+1, NY-NPOINTS_PML
      do i = NPOINTS_PML+1, NX-NPOINTS_PML
! interpolate material parameters at the right location in the staggered grid cell
        lambda_half_x = 0.5d0 * (lambda_unrelaxed(i+1,j) + lambda_unrelaxed(i,j))
        mu_half_x = 0.5d0 * (mu_unrelaxed(i+1,j) + mu_unrelaxed(i,j))
        mu_half_y = 0.5d0 * (mu_unrelaxed(i,j+1) + mu_unrelaxed(i,j))
        epsilon_xx = ((lambda_half_x + 2.d0*mu_half_x) * sigma_xx(i,j) - lambda_half_x * &
          sigma_yy(i,j)) / (4.d0 * mu_half_x * (lambda_half_x + mu_half_x))
        epsilon_yy = ((lambda_half_x + 2.d0*mu_half_x) * sigma_yy(i,j) - lambda_half_x * &
          sigma_xx(i,j)) / (4.d0 * mu_half_x * (lambda_half_x + mu_half_x))
        epsilon_xy = sigma_xy(i,j) / (2.d0 * mu_half_y)
        total_energy_potential(it) = total_energy_potential(it) + &
          0.5d0 * (epsilon_xx * sigma_xx(i,j) + epsilon_yy * sigma_yy(i,j) + 2.d0 * epsilon_xy * sigma_xy(i,j))
      enddo
    enddo

  endif

! output information
  if (mod(it,IT_DISPLAY) == 0 .or. it == 5) then

! print maximum of norm of velocity
    velocnorm = maxval(sqrt(vx**2 + vy**2))
    print *,'Time step # ',it,' out of ',NSTEP
    print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
    print *,'Max norm velocity vector V (m/s) = ',velocnorm
    if (COMPUTE_ENERGY) print *,'total energy = ',total_energy_kinetic(it) + total_energy_potential(it)
    print *
! check stability of the code, exit if unstable
    if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

    call create_color_image(vx,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1)
    call create_color_image(vy,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2)

! save the part of the seismograms that has been computed so far, so that users can monitor the progress of the simulation
    call write_seismograms(sisvx,sisvy,sispressure,NSTEP,NREC,DELTAT,t0)

  endif

  enddo   ! end of time loop

! save seismograms
  call write_seismograms(sisvx,sisvy,sispressure,NSTEP,NREC,DELTAT,t0)

  if (COMPUTE_ENERGY) then

! save total energy
    open(unit=20,file='energy.dat',status='unknown')
    do it = 1,NSTEP
      write(20,*) sngl(dble(it-1)*DELTAT),sngl(total_energy_kinetic(it)), &
         sngl(total_energy_potential(it)),sngl(total_energy_kinetic(it) + total_energy_potential(it))
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
    write(20,*) 'set output "cpml_total_energy_semilog.eps"'
    write(20,*) 'set logscale y'
    write(20,*) 'plot "energy.dat" us 1:2 t ''Ec'' w l lc 1, "energy.dat" us 1:3 &
                & t ''Ep'' w l lc 3, "energy.dat" us 1:4 t ''Total energy'' w l lc 4'
    write(20,*) 'pause -1 "Hit any key..."'
    write(20,*)
    close(20)

  endif

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

  write(20,*) 'set output "v_sigma_Vx_receiver_002.eps"'
  write(20,*) 'plot "Vx_file_002.dat" t ''Vx C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_002.eps"'
  write(20,*) 'plot "Vy_file_002.dat" t ''Vy C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  close(20)

  print *
  print *,'End of the simulation'
  print *

  end program seismic_CPML_2D_viscoelast_fourth

!----
!----  save the seismograms in ASCII text format
!----

  subroutine write_seismograms(sisvx,sisvy,sispressure,nt,nrec,DELTAT,t0)

  implicit none

  integer nt,nrec
  double precision DELTAT,t0

  double precision sisvx(nt,nrec)
  double precision sisvy(nt,nrec)
  double precision sispressure(nt,nrec)

  integer irec,it

  character(len=100) file_name

! pressure
  do irec=1,nrec
    write(file_name,"('pressure_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
! in the scheme of eq (13) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
! pressure is defined at time t + DELTAT/2, i.e. staggered in time with respect to velocity.
! Here we must thus take this shift of DELTAT/2 into account to save the seismograms at the right time
      write(11,*) sngl(dble(it-1)*DELTAT - t0 + DELTAT/2.d0),' ',sngl(sispressure(it,irec))
    enddo
    close(11)
  enddo

! X component of velocity
  do irec=1,nrec
    write(file_name,"('Vx_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvx(it,irec))
    enddo
    close(11)
  enddo

! Y component of velocity
  do irec=1,nrec
    write(file_name,"('Vy_file_half_a_grid_cell_away_from_Vx_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvy(it,irec))
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

! in order to be able to use a fourth-order spatial operator on the edges of the model
! here we define the array with size (0:NX+1,0:NY+1) instead of size (NX,NY) as in the second-order case
  double precision, dimension(0:NX+1,0:NY+1) :: image_data_2D

  integer, dimension(nrec) :: ix_rec,iy_rec

  integer :: ix,iy,irec

  character(len=100) :: file_name,system_command

  integer :: R, G, B

  double precision :: normalized_value,max_amplitude

! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org
  if (field_number == 1) then
    write(file_name,"('image',i6.6,'_Vx.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vx.pnm image',i6.6,'_Vx.gif ; rm image',i6.6,'_Vx.pnm')") it,it,it
  else if (field_number == 2) then
    write(file_name,"('image',i6.6,'_Vy.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vy.pnm image',i6.6,'_Vy.gif ; rm image',i6.6,'_Vy.pnm')") it,it,it
  else if (field_number == 3) then
    write(file_name,"('image',i6.6,'_pressure.pnm')") it
    write(system_command,"('convert image',i6.6,'_pressure.pnm image',i6.6,'_pressure.gif ; rm image',i6.6,'_pressure.pnm')") &
                               it,it,it
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

!
!---- include the SolvOpt() routine that is used to compute the tau_epsilon and tau_sigma values from a given Q attenuation factor
!

include "attenuation_model_with_SolvOpt.f90"

