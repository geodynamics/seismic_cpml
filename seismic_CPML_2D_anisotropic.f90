!
! SEISMIC_CPML Version 1.1.1, November 2009.
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, komatitsch aT lma DOT cnrs-mrs DOT fr
!               and Roland Martin, roland DOT martin aT get DOT obs-mip DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional anisotropic elastic wave equation
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

  program seismic_CPML_2D_aniso

! 2D elastic finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions for an anisotropic medium

! Dimitri Komatitsch, University of Pau, France, April 2007.
! Anisotropic implementation by Roland Martin and Dimitri Komatitsch, University of Pau, France, April 2007.

! The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
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

! The C-PML implementation is based in part on formulas given in Roden and Gedney (2000).
! If you use this code for your own research, please cite some (or all) of these
! articles:
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
! If you use the anisotropic implementation, please cite this article,
! in which the anisotropic parameters are described, as well:
!
! @ARTICLE{KoBaTr00,
! author = {D. Komatitsch and C. Barnes and J. Tromp},
! title = {Simulation of anisotropic wave propagation based upon a spectral element method},
! journal = {Geophysics},
! year = {2000},
! volume = {65},
! number = {4},
! pages = {1251-1260},
! doi = {10.1190/1.1444816}}
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

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 401
  integer, parameter :: NY = 401

! size of a grid cell
  double precision, parameter :: DELTAX = 0.0625d-2
  double precision, parameter :: DELTAY = DELTAX

! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 10

! Velocity of qP along horizontal axis  = sqrt(c11/rho)
! Velocity of qP along vertical axis    = sqrt(c22/rho)
! Velocity of qSV along horizontal axis = sqrt(c33/rho)
! Velocity of qSV along vertical axis   = sqrt(c33/rho), same as along horizontal axis

! zinc, from Komatitsch et al. (2000)
! double precision, parameter :: c11 = 16.5d10
! double precision, parameter :: c12 = 5.d10
! double precision, parameter :: c22 = 6.2d10
! double precision, parameter :: c33 = 3.96d10
! double precision, parameter :: rho = 7100.d0
! double precision, parameter :: f0 = 170.d3

! apatite, from Komatitsch et al. (2000)
! double precision, parameter :: c11 = 16.7d10
! double precision, parameter :: c12 = 6.6d10
! double precision, parameter :: c22 = 14.d10
! double precision, parameter :: c33 = 6.63d10
! double precision, parameter :: rho = 3200.d0
! double precision, parameter :: f0 = 300.d3

! isotropic material a bit similar to apatite
! double precision, parameter :: c11 = 16.7d10
! double precision, parameter :: c12 = c11/3.d0
! double precision, parameter :: c22 = c11
! double precision, parameter :: c33 = (c11-c12)/2.d0  ! = c11/3.d0
! double precision, parameter :: rho = 3200.d0
! double precision, parameter :: f0 = 300.d3

! model I from Becache, Fauqueux and Joly, which is stable
  double precision, parameter :: scale_aniso = 1.d10
  double precision, parameter :: c11 = 4.d0 * scale_aniso
  double precision, parameter :: c12 = 3.8d0 * scale_aniso
  double precision, parameter :: c22 = 20.d0 * scale_aniso
  double precision, parameter :: c33 = 2.d0 * scale_aniso
  double precision, parameter :: rho = 4000.d0  ! used to be 1.
  double precision, parameter :: f0 = 200.d3

! model II from Becache, Fauqueux and Joly, which is stable
! double precision, parameter :: scale_aniso = 1.d10
! double precision, parameter :: c11 = 20.d0 * scale_aniso
! double precision, parameter :: c12 = 3.8d0 * scale_aniso
! double precision, parameter :: c22 = c11
! double precision, parameter :: c33 = 2.d0 * scale_aniso
! double precision, parameter :: rho = 4000.d0  ! used to be 1.
! double precision, parameter :: f0 = 200.d3

! model III from Becache, Fauqueux and Joly, which is unstable
! double precision, parameter :: scale_aniso = 1.d10
! double precision, parameter :: c11 = 4.d0 * scale_aniso
! double precision, parameter :: c12 = 4.9d0 * scale_aniso
! double precision, parameter :: c22 = 20.d0 * scale_aniso
! double precision, parameter :: c33 = 2.d0 * scale_aniso
! double precision, parameter :: rho = 4000.d0  ! used to be 1.
! double precision, parameter :: f0 = 250.d3

! model IV from Becache, Fauqueux and Joly, which is unstable
! double precision, parameter :: scale_aniso = 1.d10
! double precision, parameter :: c11 = 4.d0 * scale_aniso
! double precision, parameter :: c12 = 7.5d0 * scale_aniso
! double precision, parameter :: c22 = 20.d0 * scale_aniso
! double precision, parameter :: c33 = 2.d0 * scale_aniso
! double precision, parameter :: rho = 4000.d0  ! used to be 1.
! double precision, parameter :: f0 = 170.d3

! total number of time steps
  integer, parameter :: NSTEP = 3000

! time step in seconds
  double precision, parameter :: DELTAT = 50.d-9

! parameters for the source
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.d7

! source
  integer, parameter :: ISOURCE = NX / 2
  integer, parameter :: JSOURCE = NY / 2
  double precision, parameter :: xsource = (ISOURCE - 1) * DELTAX
  double precision, parameter :: ysource = (JSOURCE - 1) * DELTAY
! angle of source force clockwise with respect to vertical (Y) axis
  double precision, parameter :: ANGLE_FORCE = 0.d0

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
  double precision, dimension(NX,NY) :: vx,vy,sigmaxx,sigmayy,sigmaxy

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
      memory_dsigmaxx_dx, &
      memory_dsigmayy_dy, &
      memory_dsigmaxy_dx, &
      memory_dsigmaxy_dy

  double precision :: &
      value_dvx_dx, &
      value_dvx_dy, &
      value_dvy_dx, &
      value_dvy_dy, &
      value_dsigmaxx_dx, &
      value_dsigmayy_dy, &
      value_dsigmaxy_dx, &
      value_dsigmaxy_dy

! 1D arrays for the damping profiles
  double precision, dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half
  double precision, dimension(NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half

  double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
  double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

! for the source
  double precision :: a,t,force_x,force_y,source_term

  integer :: i,j,it

  double precision :: Courant_number,velocnorm

! for stability estimate
  double precision :: quasi_cp_max,aniso_stability_criterion,aniso2,aniso3

!---
!--- program starts here
!---

  print *
  print *,'2D elastic finite-difference code in velocity and stress formulation with C-PML'
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

  print *,'Velocity of qP along vertical axis. . . . =',sqrt(c22/rho)
  print *,'Velocity of qP along horizontal axis. . . =',sqrt(c11/rho)
  print *
  print *,'Velocity of qSV along vertical axis . . . =',sqrt(c33/rho)
  print *,'Velocity of qSV along horizontal axis . . =',sqrt(c33/rho)
  print *

! from Becache et al., INRIA report, equation 7 page 5 http://hal.inria.fr/docs/00/07/22/83/PDF/RR-4304.pdf
  if (c11*c22 - c12*c12 <= 0.d0) stop 'problem in definition of orthotropic material'

! check intrinsic mathematical stability of PML model for an anisotropic material
! from E. B\'ecache, S. Fauqueux and P. Joly, Stability of Perfectly Matched Layers, group
! velocities and anisotropic waves, Journal of Computational Physics, 188(2), p. 399-433 (2003)
  aniso_stability_criterion = ((c12+c33)**2 - c11*(c22-c33)) * ((c12+c33)**2 + c33*(c22-c33))
  print *,'PML anisotropy stability criterion from Becache et al. 2003 = ',aniso_stability_criterion
  if (aniso_stability_criterion > 0.d0 .and. (USE_PML_XMIN .or. USE_PML_XMAX .or. USE_PML_YMIN .or. USE_PML_YMAX)) &
     print *,'WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 1'
  print *

  aniso2 = (c12 + 2*c33)**2 - c11*c22
  print *,'PML aniso2 stability criterion from Becache et al. 2003 = ',aniso2
  if (aniso2 > 0.d0 .and. (USE_PML_XMIN .or. USE_PML_XMAX .or. USE_PML_YMIN .or. USE_PML_YMAX)) &
     print *,'WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 2'
  print *

  aniso3 = (c12 + c33)**2 - c11*c22 - c33**2
  print *,'PML aniso3 stability criterion from Becache et al. 2003 = ',aniso3
  if (aniso3 > 0.d0 .and. (USE_PML_XMIN .or. USE_PML_XMAX .or. USE_PML_YMIN .or. USE_PML_YMAX)) &
     print *,'WARNING: PML model mathematically intrinsically unstable for this anisotropic material for condition 3'
  print *

! to compute d0 below, and for stability estimate
  quasi_cp_max = max(sqrt(c22/rho),sqrt(c11/rho))

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
  thickness_PML_x = NPOINTS_PML * DELTAX
  thickness_PML_y = NPOINTS_PML * DELTAY

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0

! check that NPOWER is okay
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * quasi_cp_max * log(Rcoef) / (2.d0 * thickness_PML_y)

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

! print position of the source
  print *,'Position of the source:'
  print *
  print *,'x = ',xsource
  print *,'y = ',ysource
  print *

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
  Courant_number = quasi_cp_max * DELTAT * sqrt(1.d0/DELTAX**2 + 1.d0/DELTAY**2)
  print *,'Courant number is ',Courant_number
  print *
  if (Courant_number > 1.d0) stop 'time step is too large, simulation will be unstable'

! suppress old files (can be commented out if "call system" is missing in your compiler)
! call system('rm -f Vx_*.dat Vy_*.dat image*.pnm image*.gif')

! initialize arrays
  vx(:,:) = ZERO
  vy(:,:) = ZERO
  sigmaxx(:,:) = ZERO
  sigmayy(:,:) = ZERO
  sigmaxy(:,:) = ZERO

! PML
  memory_dvx_dx(:,:) = ZERO
  memory_dvx_dy(:,:) = ZERO
  memory_dvy_dx(:,:) = ZERO
  memory_dvy_dy(:,:) = ZERO
  memory_dsigmaxx_dx(:,:) = ZERO
  memory_dsigmayy_dy(:,:) = ZERO
  memory_dsigmaxy_dx(:,:) = ZERO
  memory_dsigmaxy_dy(:,:) = ZERO

!---
!---  beginning of time loop
!---

  do it = 1,NSTEP

!------------------------------------------------------------
! compute stress sigma and update memory variables for C-PML
!------------------------------------------------------------

  do j = 2,NY
    do i = 1,NX-1

      value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DELTAX
      value_dvy_dy = (vy(i,j) - vy(i,j-1)) / DELTAY

      memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx
      memory_dvy_dy(i,j) = b_y(j) * memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy

      value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
      value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j)

      sigmaxx(i,j) = sigmaxx(i,j) + (c11 * value_dvx_dx + c12 * value_dvy_dy) * DELTAT
      sigmayy(i,j) = sigmayy(i,j) + (c12 * value_dvx_dx + c22 * value_dvy_dy) * DELTAT

    enddo
  enddo

  do j = 1,NY-1
    do i = 2,NX

      value_dvy_dx = (vy(i,j) - vy(i-1,j)) / DELTAX
      value_dvx_dy = (vx(i,j+1) - vx(i,j)) / DELTAY

      memory_dvy_dx(i,j) = b_x(i) * memory_dvy_dx(i,j) + a_x(i) * value_dvy_dx
      memory_dvx_dy(i,j) = b_y_half(j) * memory_dvx_dy(i,j) + a_y_half(j) * value_dvx_dy

      value_dvy_dx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j)
      value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j)

      sigmaxy(i,j) = sigmaxy(i,j) + c33 * (value_dvy_dx + value_dvx_dy) * DELTAT

    enddo
  enddo

!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

  do j = 2,NY
    do i = 2,NX

      value_dsigmaxx_dx = (sigmaxx(i,j) - sigmaxx(i-1,j)) / DELTAX
      value_dsigmaxy_dy = (sigmaxy(i,j) - sigmaxy(i,j-1)) / DELTAY

      memory_dsigmaxx_dx(i,j) = b_x(i) * memory_dsigmaxx_dx(i,j) + a_x(i) * value_dsigmaxx_dx
      memory_dsigmaxy_dy(i,j) = b_y(j) * memory_dsigmaxy_dy(i,j) + a_y(j) * value_dsigmaxy_dy

      value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j)
      value_dsigmaxy_dy = value_dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j)

      vx(i,j) = vx(i,j) + (value_dsigmaxx_dx + value_dsigmaxy_dy) * DELTAT / rho

    enddo
  enddo

  do j = 1,NY-1
    do i = 1,NX-1

      value_dsigmaxy_dx = (sigmaxy(i+1,j) - sigmaxy(i,j)) / DELTAX
      value_dsigmayy_dy = (sigmayy(i,j+1) - sigmayy(i,j)) / DELTAY

      memory_dsigmaxy_dx(i,j) = b_x_half(i) * memory_dsigmaxy_dx(i,j) + a_x_half(i) * value_dsigmaxy_dx
      memory_dsigmayy_dy(i,j) = b_y_half(j) * memory_dsigmayy_dy(i,j) + a_y_half(j) * value_dsigmayy_dy

      value_dsigmaxy_dx = value_dsigmaxy_dx / K_x_half(i) + memory_dsigmaxy_dx(i,j)
      value_dsigmayy_dy = value_dsigmayy_dy / K_y_half(j) + memory_dsigmayy_dy(i,j)

      vy(i,j) = vy(i,j) + (value_dsigmaxy_dx + value_dsigmayy_dy) * DELTAT / rho

    enddo
  enddo

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

  vx(i,j) = vx(i,j) + force_x * DELTAT / rho
  vy(i,j) = vy(i,j) + force_y * DELTAT / rho

! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  vx(1,:) = ZERO
  vx(NX,:) = ZERO

  vx(:,1) = ZERO
  vx(:,NY) = ZERO

  vy(1,:) = ZERO
  vy(NX,:) = ZERO

  vy(:,1) = ZERO
  vy(:,NY) = ZERO

! output information
  if (mod(it,IT_DISPLAY) == 0 .or. it == 5) then

! print maximum of norm of velocity
    velocnorm = maxval(sqrt(vx**2 + vy**2))
    print *,'Time step # ',it,' out of ',NSTEP
    print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
    print *,'Max norm velocity vector V (m/s) = ',velocnorm
    print *
! check stability of the code, exit if unstable
    if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

    call create_color_image(vx,NX,NY,it,ISOURCE,JSOURCE, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1)
    call create_color_image(vy,NX,NY,it,ISOURCE,JSOURCE, &
                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2)

  endif

  enddo   ! end of time loop

  print *
  print *,'End of the simulation'
  print *

  end program seismic_CPML_2D_aniso

!----
!----  routine to create a color image of a given vector component
!----  the image is created in PNM format and then converted to GIF
!----

  subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE, &
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

  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX

  double precision, dimension(NX,NY) :: image_data_2D

  integer :: ix,iy

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

! write color pixel
    write(27,"(i3,' ',i3,' ',i3)") R,G,B

    enddo
  enddo

! close file
  close(27)

! call the system to convert image to Gif (can be commented out if "call system" is missing in your compiler)
! call system(system_command)

  end subroutine create_color_image

