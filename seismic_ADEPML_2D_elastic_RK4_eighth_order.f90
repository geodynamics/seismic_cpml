!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Roland Martin, roland DOT martin aT get DOT obs-mip DOT fr
!               and Youshan Liu, China.
!
! This software is a computer program whose purpose is to solve
! the two-dimensional isotropic elastic wave equation
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

program seismic_ADEPML_2D_elastic_RK4_eighth_order

! High order 2D explicit-semi implicit-implicit elastic finite-difference code
! in velocity and stress formulation with Auxiliary Differential
! Equation Perfectly Matched Layer (ADE-PML) absorbing conditions for
! an isotropic elastic medium. It is fourth order Runge-Kutta (RK4) in time
! and 8th order in space using Holberg spatial discretization.

! Version 1.1.3
! by Roland Martin, University of Pau, France, Jan 2010
! with a major bug fix in the Runge-Kutta implementation
! and also significant memory usage optimization by Youshan Liu, China, August 2015.
! based on seismic_CPML_2D_isotropic_second_order.f90
! by Dimitri Komatitsch and Roland Martin, University of Pau, France, 2007.

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
!integer, parameter :: NX = 101
!integer, parameter :: NY = 641
integer, parameter :: NX = 241
integer, parameter :: NY = 241

! Explicit (epsn=1,epsn=0), implicit (epsn=0,epsn1=1), semi-implicit (epsn=0.5,epsn1=0.5)
integer, parameter :: iexpl=0
integer, parameter :: iimpl=0
integer, parameter :: isemiimpl=1

double precision :: epsn,epsn1

! size of a grid cell
double precision, parameter :: DELTAX = 10.d0
double precision, parameter :: DELTAY = DELTAX

! flags to add PML layers to the edges of the grid
logical, parameter :: USE_PML_XMIN = .true.
logical, parameter :: USE_PML_XMAX = .true.
logical, parameter :: USE_PML_YMIN = .true.
logical, parameter :: USE_PML_YMAX = .true.

! thickness of the PML layer in grid points. 8th order in space imposes to
! increase the thickness of the CPML
integer, parameter :: NPOINTS_PML = 10

! P-velocity, S-velocity and density
double precision, parameter :: cp      =  2000.d0
double precision, parameter :: cs      =  1150.d0
double precision, parameter :: density = 2000.d0
!double precision, parameter :: cp = 3300.d0
!double precision, parameter :: cs =  1905.d0
!double precision, parameter :: density = 2800.d0

! total number of time steps
! the time step is twice smaller for this fourth-order simulation,
! therefore let us double the number of time steps to keep the same total duration
integer, parameter :: NSTEP =  2501

! time step in seconds
! 8th-order in space and 4th-order in time finite-difference schemes
! are less stable than second-order in space and second-order in time,
! therefore let us divide the time step by 2
double precision, parameter :: DELTAT = 3.d-3

! parameters for the source
double precision, parameter :: f0 = 10.d0
double precision, parameter :: t0 = 1.0d0 / f0
double precision, parameter :: factor = 1.d4

! source
!integer, parameter :: ISOURCE = NX - 2*NPOINTS_PML  - 1
integer, parameter :: ISOURCE = (NX-1)/2
integer, parameter :: JSOURCE = (NY-1)/2
double precision, parameter :: xsource = (ISOURCE - 1) * DELTAX
double precision, parameter :: ysource = (JSOURCE - 1) * DELTAY
! angle of source force clockwise with respect to vertical (Y) axis
!double precision, parameter :: ANGLE_FORCE = 135.d0
double precision, parameter :: ANGLE_FORCE = 90.d0

! receivers
!integer, parameter :: NREC = 3
!double precision, parameter :: xdeb = xsource    ! first receiver x in meters
!double precision, parameter :: ydeb = ysource - 2000.d0   ! first receiver y in meters
!double precision, parameter :: xfin = xsource    ! last receiver x in meters
!double precision, parameter :: yfin = ysource - 4000.d0  ! last receiver y in meters
integer, parameter :: NREC = NX
double precision, parameter :: xdeb = 0.d0    ! first receiver x in meters
double precision, parameter :: ydeb = 50.d0   ! first receiver y in meters
double precision, parameter :: xfin = (NX-1)*DELTAX    ! last receiver x in meters
double precision, parameter :: yfin = 50.d0  ! last receiver y in meters

! display information on the screen from time to time
! the time step is twice smaller for this fourth-order simulation,
! therefore let us double the interval in time steps at which we display information
integer, parameter :: IT_DISPLAY = 200
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
!  @ARTICLE{Hol87,
!  author = {O. Holberg},
!  title = {Computational aspects of the choice of operator and sampling interval
!  for numerical differentiation in large-scale simulation of wave phenomena},
!  journal = {Geophysical Prospecting},
!  year = {1987},
!  volume = {35},
!  pages = {629-655}}
double precision, parameter :: c1 = 1.231666d0
double precision, parameter :: c2 = -1.041182d-1
double precision, parameter :: c3 = 2.063707d-2
double precision, parameter :: c4 = -3.570998d-3

! RK4 scheme coefficients, 2 per subloop, 8 in total
double precision, dimension(4) :: rk41, rk42

! main arrays
double precision, dimension(-4:NX+4,-4:NY+4) :: lambda,mu,rho,vx,vy,sigmaxx,sigmayy,sigmaxy

! variables are stored in four indices in the first dimension to implement RK4
! dv does not always indicate a derivative
double precision, dimension(3,-4:NX+4,-4:NY+4) :: dvx,dvy,dsigmaxx,dsigmayy,dsigmaxy

! to interpolate material parameters at the right location in the staggered grid cell
double precision lambda_half_x,mu_half_x,lambda_plus_two_mu_half_x,mu_half_y,rho_half_x_half_y

! for evolution of total energy in the medium
double precision, dimension(NSTEP) :: total_energy_kinetic,total_energy_potential

! power to compute d0 profile
double precision, parameter :: NPOWER = 2.d0
double precision, parameter :: NPOWER2 = 2.d0

! Kappa must be strong enough to absorb energy and low enough to avoid
! reflections.
! Alpha1 must be low to absorb energy and high enough to have efficiency on
! grazing incident waves.
double precision, parameter :: K_MAX_PML = 7.d0
double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0)

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
!!! Youshan Liu suppressed the two comment lines below
!!!!!! not true anymore: We have as many memory variables as the number of frequency shift poles in the CPML
!!!!!! not true anymore: Indices are 1 and 2 for the 2 frequency shift poles
! ==================== revised by Youshan Liu ==================
double precision, dimension(-4:NX+4,-4:NY+4) :: memory_dvx_dx, memory_dvx_dy, memory_dvy_dx, memory_dvy_dy, &
                                                  memory_dsigmaxx_dx, memory_dsigmayy_dy, &
                                                  memory_dsigmaxy_dx, memory_dsigmaxy_dy

double precision :: value_dvx_dx, value_dvx_dy, value_dvy_dx, value_dvy_dy, &
                    value_dsigmaxx_dx, value_dsigmayy_dy, &
                    value_dsigmaxy_dx, value_dsigmaxy_dy

! 1D arrays for the damping profiles
double precision, dimension(-4:NX+4) :: d_x,K_x,alpha_x,g_x,ksi_x
double precision, dimension(-4:NX+4) :: d_x_half,K_x_half,alpha_x_half,g_x_half,ksi_x_half
double precision, dimension(-4:NY+4) :: d_y,K_y,alpha_y,g_y,ksi_y
double precision, dimension(-4:NY+4) :: d_y_half,K_y_half,alpha_y_half,g_y_half,ksi_y_half

! coefficients that allow to reset the memory variables at each RK4 substep depend on the substepping and are then of dimension 4,
! 1D arrays for the damping profiles
double precision, dimension(4,-4:NX+4) :: a_x,b_x
double precision, dimension(4,-4:NX+4) :: a_x_half,b_x_half
double precision, dimension(4,-4:NY+4) :: a_y,b_y
double precision, dimension(4,-4:NY+4) :: a_y_half,b_y_half


double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

! for the source
double precision :: a,t,force_x,force_y,source_term

! for receivers
double precision xspacerec,yspacerec,distval,dist
integer, dimension(NREC) :: ix_rec,iy_rec
double precision, dimension(NREC) :: xrec,yrec

! for seismograms
double precision, dimension(NSTEP,NREC) :: sisvx,sisvy

integer :: i,j,k,it,irec,inc

double precision :: Courant_number

!define by ysliu 8/2/2015
integer(2) head(1:120)

character(80) :: routine
real,dimension(NSTEP,NREC) :: seisvx, seisvy
real,dimension(NX,NY) :: snapvx,snapvy

!---
!--- program starts here
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

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
thickness_PML_x = NPOINTS_PML * DELTAX
thickness_PML_y = NPOINTS_PML * DELTAY

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
Rcoef = 0.00001d0

! check that NPOWER is okay
if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
d0_x = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_x)
d0_y = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_y)

print *,'d0_x = ',d0_x
print *,'d0_y = ',d0_y
print *

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
d_x(:) = ZERO
d_x_half(:) = ZERO
K_x(:) = 1.d0
K_x_half(:) = 1.d0
alpha_x(:) = ZERO
alpha_x_half(:) = ZERO
a_x(:,:) = ZERO
a_x_half(:,:) = ZERO
g_x(:) = 5.d-1
g_x_half(:) = 5.d-1

ksi_y(:) = ZERO
ksi_y_half(:) = ZERO
d_y(:) = ZERO
d_y_half(:) = ZERO
K_y(:) = 1.d0
K_y_half(:) = 1.d0
alpha_y(:) = ZERO
alpha_y_half(:) = ZERO
a_y(:,:) = ZERO
a_y_half(:,:) = ZERO
g_y(:) = 1.d0
g_y_half(:) = 1.d0

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
         d_x(i) = d0_x * abscissa_normalized**NPOWER
         ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
         K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
         alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

      ! define damping profile at half the grid points
      abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0)
      if (abscissa_in_PML >= ZERO) then
         abscissa_normalized = abscissa_in_PML / thickness_PML_x
         d_x_half(i) = d0_x * abscissa_normalized**NPOWER
         ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
         K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
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
         K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
         alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

      ! define damping profile at half the grid points
      abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright
      if (abscissa_in_PML >= ZERO) then
         abscissa_normalized = abscissa_in_PML / thickness_PML_x
         d_x_half(i) = d0_x * abscissa_normalized**NPOWER
         ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
         K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
         alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

   endif

   ! just in case, for -5 at the end
   if (alpha_x(i) < ZERO) alpha_x(i) = ZERO
   if (alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

   ! CPML damping parameters for the 4 sub time steps of RK4 algorithm
   do inc=1,4
      b_x(inc,i) =  (1.-epsn*DELTAT*rk41(inc)*(d_x(i)/K_x(i) + alpha_x(i)))/ &
                         (1.+epsn1*DELTAT*rk41(inc)*(d_x(i)/K_x(i) + alpha_x(i)))
      b_x_half(inc,i) = (1.-epsn*DELTAT*rk41(inc)*(d_x_half(i)/K_x_half(i) &
           + alpha_x_half(i)))/(1. +epsn1*DELTAT*rk41(inc)*(d_x_half(i)/K_x_half(i) &
           + alpha_x_half(i)))

      ! this to avoid division by zero outside the PML
      if (abs(d_x(i)) > 1.d-6) a_x(inc,i) = - DELTAT*rk41(inc)*d_x(i) / (K_x(i)* K_x(i))/&
          (1. +epsn1*DELTAT*rk41(inc)*(d_x(i)/K_x(i) + alpha_x(i)))

      if (abs(d_x_half(i)) > 1.d-6) a_x_half(inc,i) =-DELTAT*rk41(inc)*d_x_half(i)/&
          (K_x_half(i)*K_x_half(i) )/&
          (1. +epsn1*DELTAT*rk41(inc)*(d_x_half(i)/K_x_half(i)&
          + alpha_x_half(i)))

    enddo

enddo !do i = -4,NX+4

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
         d_y(j) = d0_y * abscissa_normalized**NPOWER
         ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
         K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
         alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

      ! define damping profile at half the grid points
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0)
      if (abscissa_in_PML >= ZERO) then
         abscissa_normalized = abscissa_in_PML / thickness_PML_y
         d_y_half(j) = d0_y * abscissa_normalized**NPOWER
         ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
         K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
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
         K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
         alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

      ! define damping profile at half the grid points
      abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop
      if (abscissa_in_PML >= ZERO) then
         abscissa_normalized = abscissa_in_PML / thickness_PML_y
         d_y_half(j) = d0_y * abscissa_normalized**NPOWER
         ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
         K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER2
         alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

   endif

   ! just in case, for -5 at the end
   if (alpha_y(j) < ZERO) alpha_y(j) = ZERO
   if (alpha_y_half(j) < ZERO) alpha_y_half(j) = ZERO

   ! CPML damping parameters for the 4 sub time steps of RK4 algorithm
   do inc=1,4
      b_y(inc,j) =  (1.-epsn*DELTAT*rk41(inc)*(d_y(j)/K_y(j) + alpha_y(j)))/ &
                       (1.+epsn1*DELTAT*rk41(inc)*(d_y(j)/K_y(j) + alpha_y(j)))
      b_y_half(inc,j) = (1.-epsn*DELTAT*rk41(inc)*(d_y_half(j)/K_y_half(j) + &
      alpha_y_half(j)))/(1.+epsn1*DELTAT*rk41(inc)*(d_y_half(j)/K_y_half(j)  &
                                                                + alpha_y_half(j)))

      ! this to avoid division by zero outside the PML
      if (abs(d_y(j)) > 1.d-6) a_y(inc,j) = - DELTAT*rk41(inc)*d_y(j) &
                  / (K_y(j)* K_y(j))/&
                  (1.+epsn1*DELTAT*rk41(inc)*(d_y(j)/K_y(j) + alpha_y(j)))
      if (abs(d_y_half(j)) > 1.d-6) a_y_half(inc,j) = -DELTAT*rk41(inc)*d_y_half(j) /&
                  (K_y_half(j) * K_y_half(j) )/&
                  (1.+epsn1*DELTAT*rk41(inc)*(d_y_half(j)/K_y_half(j) + alpha_y_half(j)))
   enddo

enddo !do j = -4,NY+4

! compute the Lame parameters and density
do j = -4,NY+4
   do i = -4,NX+4
      rho(i,j) = density
      mu(i,j) = density*cs*cs
      lambda(i,j) = density*(cp*cp - 2.d0*cs*cs)
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
xspacerec = (xfin-xdeb) / dble(NREC-1)
yspacerec = (yfin-ydeb) / dble(NREC-1)
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
enddo !do irec=1,nrec

! check the Courant stability condition for the explicit time scheme
! R. Courant and K. O. Friedrichs and H. Lewy (1928)
Courant_number = cp * DELTAT * sqrt(1.d0/DELTAX**2 + 1.d0/DELTAY**2)
print *,'Courant number is ',Courant_number
print *
if (Courant_number > 1.d0) stop 'time step is too large, simulation will be unstable'

! suppress old files (can be commented out if "call system" is missing in your compiler)
! call system('rm -f Vx_*.dat Vy_*.dat image*.pnm image*.gif')

! initialize arrays
dvx(:,:,:) = ZERO
dvy(:,:,:) = ZERO
dsigmaxx(:,:,:) = ZERO
dsigmayy(:,:,:) = ZERO
dsigmaxy(:,:,:) = ZERO

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

! initialize seismograms
sisvx(:,:) = ZERO
sisvy(:,:) = ZERO

! initialize total energy
total_energy_kinetic(:) = ZERO
total_energy_potential(:) = ZERO

!---
!---  beginning of time loop
!---

do it = 1,NSTEP

   !! v and sigma temporary variables of RK4
   !======================================================
   !====================revised by ysliu==================
   !backup the current snapshots
   dvx(2,:,:) = vx(:,:)
   dvy(2,:,:) = vy(:,:)
   dsigmaxx(2,:,:) = sigmaxx(:,:)
   dsigmayy(2,:,:) = sigmayy(:,:)
   dsigmaxy(2,:,:) = sigmaxy(:,:)
   dvx(3,:,:) = vx(:,:)
   dvy(3,:,:) = vy(:,:)
   dsigmaxx(3,:,:) = sigmaxx(:,:)
   dsigmayy(3,:,:) = sigmayy(:,:)
   dsigmaxy(3,:,:) = sigmaxy(:,:)

   !======================================================

   ! RK4 loop (loop on the four RK4 substeps)
   do inc= 1,4
      ! ==================== revised by Youshan Liu ==================
      ! The new values of the different variables v and sigma are computed
      dvx(1,:,:) = dvx(3,:,:) + rk41(inc) * dvx(2,:,:) * DELTAT
      dvy(1,:,:) = dvy(3,:,:) + rk41(inc) * dvy(2,:,:) * DELTAT
      dsigmaxx(1,:,:) = dsigmaxx(3,:,:) + rk41(inc) * dsigmaxx(2,:,:) * DELTAT
      dsigmayy(1,:,:) = dsigmayy(3,:,:) + rk41(inc) * dsigmayy(2,:,:) * DELTAT
      dsigmaxy(1,:,:) = dsigmaxy(3,:,:) + rk41(inc) * dsigmaxy(2,:,:) * DELTAT

      !------------------
      ! compute velocity
      !------------------

      do j = 2,NY
         do i = 2,NX

            value_dsigmaxx_dx = ( c1 * (dsigmaxx(1,i,j) - dsigmaxx(1,i-1,j)) + c2 * (dsigmaxx(1,i+1,j) - dsigmaxx(1,i-2,j)) + &
               c3 * (dsigmaxx(1,i+2,j) - dsigmaxx(1,i-3,j)) + c4 * (dsigmaxx(1,i+3,j) - dsigmaxx(1,i-4,j)) )/ DELTAX

            value_dsigmaxy_dy = ( c1 * (dsigmaxy(1,i,j) - dsigmaxy(1,i,j-1)) + c2* (dsigmaxy(1,i,j+1) - dsigmaxy(1,i,j-2)) + &
               c3 * (dsigmaxy(1,i,j+2) - dsigmaxy(1,i,j-3)) + c4 * (dsigmaxy(1,i,j+3) - dsigmaxy(1,i,j-4)) )/ DELTAY

            if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then
               ! ==================== revised by Youshan Liu ==================
               memory_dsigmaxx_dx(i,j) = b_x(inc,i) * memory_dsigmaxx_dx(i,j) + a_x(inc,i) * value_dsigmaxx_dx
               memory_dsigmaxy_dy(i,j) = b_y(inc,j) * memory_dsigmaxy_dy(i,j) + a_y(inc,j) * value_dsigmaxy_dy

               value_dsigmaxx_dx = value_dsigmaxx_dx / K_x(i) + memory_dsigmaxx_dx(i,j)
               value_dsigmaxy_dy = value_dsigmaxy_dy / K_y(j) + memory_dsigmaxy_dy(i,j)
            endif

            dvx(2,i,j) = (value_dsigmaxx_dx + value_dsigmaxy_dy) / rho(i,j)

         enddo
      enddo

      do j = 1,NY-1
         do i = 1,NX-1

            ! interpolate density at the right location in the staggered grid cell
            rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))

            value_dsigmaxy_dx = ( c1 * (dsigmaxy(1,i+1,j) - dsigmaxy(1,i,j)) + c2 * (dsigmaxy(1,i+2,j) - dsigmaxy(1,i-1,j)) + &
               c3 * (dsigmaxy(1,i+3,j) - dsigmaxy(1,i-2,j)) + c4 * (dsigmaxy(1,i+4,j) - dsigmaxy(1,i-3,j)) )/ DELTAX

            value_dsigmayy_dy = ( c1 * (dsigmayy(1,i,j+1) - dsigmayy(1,i,j)) + c2 * (dsigmayy(1,i,j+2) - dsigmayy(1,i,j-1)) + &
               c3 * (dsigmayy(1,i,j+3) - dsigmayy(1,i,j-2)) + c4 * (dsigmayy(1,i,j+4) - dsigmayy(1,i,j-3)) )/ DELTAY

            if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then
         ! ==================== revised by Youshan Liu ==================
               memory_dsigmaxy_dx(i,j) = b_x_half(inc,i) * memory_dsigmaxy_dx(i,j) + a_x_half(inc,i) * value_dsigmaxy_dx
               memory_dsigmayy_dy(i,j) = b_y_half(inc,j) * memory_dsigmayy_dy(i,j) + a_y_half(inc,j) * value_dsigmayy_dy

               value_dsigmaxy_dx = value_dsigmaxy_dx/K_x_half(i)+memory_dsigmaxy_dx(i,j)
               value_dsigmayy_dy = value_dsigmayy_dy/K_y_half(j)+memory_dsigmayy_dy(i,j)
            endif

            dvy(2,i,j) = (value_dsigmaxy_dx + value_dsigmayy_dy) / rho_half_x_half_y

         enddo
      enddo

      ! add the source (force vector located at a given grid point)
      a = pi*pi*f0*f0
      t = (dble(it-1)+ rk41(inc)) * DELTAT

      ! Gaussian
      ! source_term = factor * exp(-a*(t-t0)**2) !

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
      rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))

      dvx(2,i,j) = dvx(2,i,j) + force_x  / rho(i,j)
      dvy(2,i,j) = dvy(2,i,j) + force_y / rho_half_x_half_y

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

      do j = 2,NY
         do i = 1,NX-1

            ! interpolate material parameters at the right location in the staggered grid cell
            lambda_half_x = 0.5d0 * (lambda(i+1,j) + lambda(i,j))
            mu_half_x = 0.5d0 * (mu(i+1,j) + mu(i,j))
            lambda_plus_two_mu_half_x = lambda_half_x + 2.d0 * mu_half_x

            value_dvx_dx = ( c1 * (dvx(1,i+1,j) - dvx(1,i,j)) + c2 * (dvx(1,i+2,j) - dvx(1,i-1,j)) + &
               c3 * (dvx(1,i+3,j) - dvx(1,i-2,j)) + c4 * (dvx(1,i+4,j) - dvx(1,i-3,j)) )/ DELTAX

            value_dvy_dy = ( c1 * (dvy(1,i,j) - dvy(1,i,j-1)) + c2 * (dvy(1,i,j+1) - dvy(1,i,j-2)) + &
               c3 * (dvy(1,i,j+2) - dvy(1,i,j-3)) + c4 * (dvy(1,i,j+3) - dvy(1,i,j-4)) )/ DELTAY

            if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then
               ! ==================== revised by Youshan Liu ==================
               memory_dvx_dx(i,j) = b_x_half(inc,i) * memory_dvx_dx(i,j) + a_x_half(inc,i) * value_dvx_dx
               memory_dvy_dy(i,j) = b_y(inc,j) * memory_dvy_dy(i,j) + a_y(inc,j) * value_dvy_dy

               value_dvx_dx = value_dvx_dx / K_x_half(i)  + memory_dvx_dx(i,j)
               value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j)
            endif

            dsigmaxx(2,i,j) = (lambda_plus_two_mu_half_x * value_dvx_dx + lambda_half_x * value_dvy_dy)
            dsigmayy(2,i,j) =  (lambda_half_x * value_dvx_dx + lambda_plus_two_mu_half_x * value_dvy_dy)

         enddo
      enddo

      do j = 1,NY-1
         do i = 2,NX

            ! interpolate material parameters at the right location in the staggered grid cell
            mu_half_y = 0.5d0 * (mu(i,j+1) + mu(i,j))

            value_dvx_dy = ( c1 * (dvx(1,i,j+1) - dvx(1,i,j)) + c2 * (dvx(1,i,j+2) - dvx(1,i,j-1)) +  &
               c3 * (dvx(1,i,j+3) - dvx(1,i,j-2)) + c4 * (dvx(1,i,j+4) - dvx(1,i,j-3)) )/ DELTAY
            value_dvy_dx = ( c1 * (dvy(1,i,j) - dvy(1,i-1,j)) + c2 * (dvy(1,i+1,j) - dvy(1,i-2,j)) + &
               c3 * (dvy(1,i+2,j) - dvy(1,i-3,j)) + c4 * (dvy(1,i+3,j) - dvy(1,i-4,j)) )/ DELTAX

            if (i <= NPOINTS_PML+2 .or. i >= NX-NPOINTS_PML-2 .or. j <= NPOINTS_PML+2 .or. j >= NY-NPOINTS_PML-2) then
               ! ==================== revised by Youshan Liu ==================
               memory_dvy_dx(i,j) = b_x(inc,i) * memory_dvy_dx(i,j) + a_x(inc,i) * value_dvy_dx
               memory_dvx_dy(i,j) = b_y_half(inc,j) * memory_dvx_dy(i,j) + a_y_half(inc,j) * value_dvx_dy

               value_dvy_dx = value_dvy_dx / K_x(i)  + memory_dvy_dx(i,j)
               value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j)
            endif

            dsigmaxy(2,i,j) = mu_half_y * (value_dvy_dx + value_dvx_dy)

          enddo
      enddo

      ! ==================== revised by Youshan Liu ==================
      ! the new values of the different variables v and sigma are computed
      vx(:,:) = vx(:,:) + rk42(inc) * dvx(2,:,:) * DELTAT
      vy(:,:) = vy(:,:) + rk42(inc) * dvy(2,:,:) * DELTAT
      sigmaxx(:,:) = sigmaxx(:,:) + rk42(inc) * dsigmaxx(2,:,:) * DELTAT
      sigmayy(:,:) = sigmayy(:,:) + rk42(inc) * dsigmayy(2,:,:) * DELTAT
      sigmaxy(:,:) = sigmaxy(:,:) + rk42(inc) * dsigmaxy(2,:,:) * DELTAT

      !! Dirichlet conditions (rigid boundaries) on all the edges of the grid
      vx(-4:1,:) = ZERO
      vx(:,-4:1) = ZERO
      vy(-4:1,:) = ZERO
      vy(:,-4:1) = ZERO

      vx(NX:NX+4,:) = ZERO
      vx(:,NY:NY+4) = ZERO
      vy(NX:NX+4,:) = ZERO
      vy(:,NY:NY+4) = ZERO

   enddo
   ! end of RK4 loop

   ! store seismograms
   do irec = 1,NREC
      sisvx(it,irec) = (vx(ix_rec(irec),iy_rec(irec))+ &
                       vx(ix_rec(irec)+1,iy_rec(irec))+ &
                       vx(ix_rec(irec),iy_rec(irec)+1)+ &
                       vx(ix_rec(irec)+1,iy_rec(irec)+1))/4.d0
      sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec))
   enddo

   !! compute total energy in the medium (without the PML layers)
   !
   !! compute kinetic energy first, defined as 1/2 rho ||v||^2
   !! in principle we should use rho_half_x_half_y instead of rho for vy
   !! in order to interpolate density at the right location in the staggered grid cell
   !! but in a homogeneous medium we can safely ignore it
   !  total_energy_kinetic(it) = 0.5d0 * sum( &
   !      rho(NPOINTS_PML+1:NX-NPOINTS_PML,NPOINTS_PML+1:NY-NPOINTS_PML)*( &
   !       vx(NPOINTS_PML+1:NX-NPOINTS_PML,NPOINTS_PML+1:NY-NPOINTS_PML)**2 +  &
   !       vy(NPOINTS_PML+1:NX-NPOINTS_PML,NPOINTS_PML+1:NY-NPOINTS_PML)**2))
   !
   !! add potential energy, defined as 1/2 epsilon_ij sigma_ij
   !! in principle we should interpolate the medium parameters at the right location
   !! in the staggered grid cell but in a homogeneous medium we can safely ignore it
   !  total_energy_potential(it) = ZERO
   !  do j = NPOINTS_PML+1, NY-NPOINTS_PML
   !    do i = NPOINTS_PML+1, NX-NPOINTS_PML
   !      epsilon_xx = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmaxx(i,j) - lambda(i,j) * &
   !        sigmayy(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)))
   !      epsilon_yy = ((lambda(i,j) + 2.d0*mu(i,j)) * sigmayy(i,j) - lambda(i,j) * &
   !        sigmaxx(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)))
   !      epsilon_xy = sigmaxy(i,j) / (2.d0 * mu(i,j))
   !      total_energy_potential(it) = total_energy_potential(it) + &
   !        0.5d0 * (epsilon_xx * sigmaxx(i,j) + epsilon_yy * sigmayy(i,j) + 2.d0 * epsilon_xy * sigmaxy(i,j))
   !    enddo
   !  enddo

   if (mod(it,IT_DISPLAY) == 0) then
      write(*,*) it, ' of ', nstep
      head=0
      head(58) = NY
      head(59) = DELTAY * 1E3
      snapvx = vx(1:NX,1:NY)
      snapvy = vy(1:NX,1:NY)
      write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapVx.su'
      open(21,file=routine,access='stream')
         do j = 1,NX,1
            write(21) head,(real(snapvx(k,j)),k=1,NY)
         enddo
      close(21)
      write(routine,'(a12,i5.5,a9)') './snapshots/',it,'snapVy.su'
      open(21,file=routine,access='stream')
         do j = 1,NX,1
            write(21) head,(real(snapvy(k,j)),k=1,NY)
         enddo
      close(21)
   endif

   !! output information
   !  if (mod(it,IT_DISPLAY) == 0 .or. it == 5) then
   !
   !! print maximum of norm of velocity
   !    velocnorm = maxval(sqrt(vx**2 + vy**2))
   !    print *,'Time step # ',it
   !    print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
   !    print *,'Max norm velocity vector V (m/s) = ',velocnorm
   !    print *,'total energy = ',total_energy_kinetic(it) + total_energy_potential(it)
   !    print *
   !! check stability of the code, exit if unstable
   !    if (velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
   !
   !    call create_color_image(vx(1:NX,1:NY),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
   !                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1)
   !    call create_color_image(vy(1:NX,1:NY),NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
   !                         NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2)
   !  open(unit=20,file='energy.dat',status='unknown')
   !  do it2 = 1,NSTEP
   !    write(20,*) sngl(dble(it2-1)*DELTAT),sngl(total_energy_kinetic(it2)), &
   !       sngl(total_energy_potential(it2)),sngl(total_energy_kinetic(it2) + total_energy_potential(it2))
   !  enddo
   !  close(20)
   !  call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT)
   !
   !  endif

enddo   ! end of time loop

! save seismograms

!save seismogram in SU format
write(*,*) NREC,nstep
seisvx = sisvx
seisvy = sisvy
head=0
head(58)=nstep
head(59)=deltat*1e6
open(21,file='./seismograms/seisVx.su',access='stream')
   do j=1,NREC,1
      write(21) head,(real(seisvx(k,j)),k=1,nstep)
   enddo
close(21)
open(21,file='./seismograms/seisVy.su',access='stream')
   do j=1,NREC,1
      write(21) head,(real(seisvy(k,j)),k=1,nstep)
   enddo
close(21)
!call write_seismograms(sisvx,sisvy,NSTEP,NREC,DELTAT)

!! save total energy
!open(unit=20,file='energy.dat',status='unknown')
!   do it = 1,NSTEP
!      write(20,*) sngl(dble(it-1)*DELTAT),sngl(total_energy_kinetic(it)), &
!         sngl(total_energy_potential(it)),sngl(total_energy_kinetic(it) + total_energy_potential(it))
!   enddo
!close(20)

!! create script for Gnuplot for total energy
!  open(unit=20,file='plot_energy',status='unknown')
!  write(20,*) '# set term x11'
!  write(20,*) 'set term postscript landscape monochrome dashed "Helvetica" 22'
!  write(20,*)
!  write(20,*) 'set xlabel "Time (s)"'
!  write(20,*) 'set ylabel "Total energy"'
!  write(20,*)
!  write(20,*) 'set output "cpml_total_energy_semilog.eps"'
!  write(20,*) 'set logscale y'
!  write(20,*) 'plot "energy.dat" us 1:2 t ''Ec'' w l lc 1, "energy.dat" us 1:3 &
!              & t ''Ep'' w l lc 3, "energy.dat" us 1:4 t ''Total energy'' w l lc 4'
!  write(20,*) 'pause -1 "Hit any key..."'
!  write(20,*)
!  close(20)
!
!  open(unit=20,file='plot_comparison',status='unknown')
!  write(20,*) '# set term x11'
!  write(20,*) 'set term postscript landscape monochrome dashed "Helvetica" 22'
!  write(20,*)
!  write(20,*) 'set xlabel "Time (s)"'
!  write(20,*) 'set ylabel "Total energy"'
!  write(20,*)
!  write(20,*) 'set output "compare_total_energy_semilog.eps"'
!  write(20,*) 'set logscale y'
!  write(20,*) 'plot "energy.dat" us 1:4 t ''Total energy CPML'' w l lc 1, &
!              & "../collino/energy.dat" us 1:4 t ''Total energy Collino'' w l lc 2'
!  write(20,*) 'pause -1 "Hit any key..."'
!  write(20,*)
!  close(20)
!
!! create script for Gnuplot
!  open(unit=20,file='plotgnu',status='unknown')
!  write(20,*) 'set term x11'
!  write(20,*) '# set term postscript landscape monochrome dashed "Helvetica" 22'
!  write(20,*)
!  write(20,*) 'set xlabel "Time (s)"'
!  write(20,*) 'set ylabel "Amplitude (m / s)"'
!  write(20,*)
!
!  write(20,*) 'set output "v_sigma_Vx_receiver_001.eps"'
!  write(20,*) 'plot "Vx_file_001.dat" t ''Vx C-PML'' w l lc 1'
!  write(20,*) 'pause -1 "Hit any key..."'
!  write(20,*)
!
!  write(20,*) 'set output "v_sigma_Vy_receiver_001.eps"'
!  write(20,*) 'plot "Vy_file_001.dat" t ''Vy C-PML'' w l lc 1'
!  write(20,*) 'pause -1 "Hit any key..."'
!  write(20,*)
!
!  write(20,*) 'set output "v_sigma_Vx_receiver_002.eps"'
!  write(20,*) 'plot "Vx_file_002.dat" t ''Vx C-PML'' w l lc 1'
!  write(20,*) 'pause -1 "Hit any key..."'
!  write(20,*)
!
!  write(20,*) 'set output "v_sigma_Vy_receiver_002.eps"'
!  write(20,*) 'plot "Vy_file_002.dat" t ''Vy C-PML'' w l lc 1'
!  write(20,*) 'pause -1 "Hit any key..."'
!  write(20,*)
!
!  close(20)

print *
print *,'End of the simulation'
print *

end program seismic_ADEPML_2D_elastic_RK4_eighth_order

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

