!
! SEISMIC_CPML Version 1.1.1, November 2009.
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Roland Martin, roland DOT martin aT get DOT obs-mip DOT fr
!           and Dimitri Komatitsch, komatitsch aT lma DOT cnrs-mrs DOT fr
!
! This software is a computer program whose purpose is to solve
! the poroelastic elastic wave equation
! using a finite-difference method with Convolutional Perfectly Matched
! Layer (C-PML) conditions and Biot model with and without viscous dissipation.
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

  program seismic_CPML_2D_poroelastic_fourth

! 2D poroelastic finite-difference code in velocity and stress formulation
! with Convolution-PML (C-PML) absorbing conditions
! with and without viscous dissipation

! Roland Martin, University of Pau, France, October 2009.
! based on the elastic code of Komatitsch and Martin, 2007.

! The fourth-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
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

! To display the results as color images in the selected 2D cut plane, use:
!
!   " display image*.gif " or " gimp image*.gif "
!
! or
!
!   " montage -geometry +0+3 -rotate 90 -tile 1x21 image*Vx*.gif allfiles_Vx.gif
!   "
!   " montage -geometry +0+3 -rotate 90 -tile 1x21 image*Vy*.gif allfiles_Vy.gif
!   "
!   then " display allfiles_Vx.gif " or " gimp allfiles_Vx.gif "
!   then " display allfiles_Vy.gif " or " gimp allfiles_Vy.gif "

! To display the 2D results as PostScript vector plots with small arrows, use:
!
!   " gs vect*.ps "
!

! IMPORTANT : all our CPML codes work fine in single precision as well (which is significantly faster).
!             If you want you can thus force automatic conversion to single precision at compile time
!             or change all the declarations and constants in the code from double precision to single.

  implicit none

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 140
  integer, parameter :: NY = 620

! size of a grid cell
  double precision, parameter :: DELTAX = 0.5D0
  double precision, parameter :: DELTAY = DELTAX

! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_LEFT   = .true.
  logical, parameter :: USE_PML_RIGHT  = .true.
  logical, parameter :: USE_PML_BOTTOM = .true.
  logical, parameter :: USE_PML_TOP    = .true.

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 10

! heterogeneous model and height of the interface
  logical, parameter :: HETEROGENEOUS_MODEL = .true.
  double precision, parameter :: INTERFACE_HEIGHT =105.D0+NPOINTS_PML*DELTAY
  integer, parameter:: JINTERFACE=INT(INTERFACE_HEIGHT/DELTAY)+1
  double precision :: co,c1,c2,vtemp

! model mud saturated with water, see article by Martin and Komatitsch
  double precision, parameter :: etaokappa_bottom=0.d0
  double precision, parameter :: rmu_bottom = 5.25D09
  double precision, parameter :: phi_bottom =0.25d0
  double precision, parameter :: a_bottom = 2.49d0
  double precision, parameter :: rhos_bottom = 2588.d0
  double precision, parameter :: rhof_bottom = 952.4d0
  double precision, parameter :: rho_bottom =2179.1d0
  double precision, parameter :: rsm_bottom =9486.d0
  double precision, parameter :: alpha_bottom=0.89d0
  double precision, parameter :: rbM_bottom =7.71d09
  double precision, parameter :: rlambdao_bottom = 6.2D08
  double precision, parameter :: rlambdac_bottom =rlambdao_bottom+alpha_bottom**2*rbM_bottom
  double precision, parameter :: ro11_b=rho_bottom+phi_bottom*rhof_bottom*(a_bottom-2.d0)
  double precision, parameter :: ro12_b=phi_bottom*rhof_bottom*(1.d0-a_bottom)
  double precision, parameter :: ro22_b=a_bottom*phi_bottom*rhof_bottom
  double precision, parameter :: lambda_b=rlambdao_bottom+rbM_bottom*(alpha_bottom-phi_bottom)**2
  double precision, parameter :: R_b=rbM_bottom*phi_bottom**2
  double precision, parameter :: ga_b=rbM_bottom*phi_bottom*(alpha_bottom-phi_bottom)
  double precision, parameter :: S_b=lambda_b+2*rmu_bottom
  double precision, parameter :: c1_b=S_b*R_b-ga_b**2
  double precision, parameter :: b1_b=-S_b*ro22_b-R_b*ro11_b+2*ga_b*ro12_b
  double precision, parameter :: a1_b=ro11_b*ro22_b-ro12_b**2
  double precision, parameter :: delta_b=b1_b**2-4.d0*a1_b*c1_b

  double precision:: cp_bottom
  double precision:: cps_bottom
  double precision:: cs_bottom

  double precision, parameter :: etaokappa_top=3.33D06
  double precision, parameter :: rmu_top = 2.4D09
  double precision, parameter :: phi_top =0.1d0
  double precision, parameter :: a_top = 2.42d0
  double precision, parameter :: rhos_top = 2250.d0
  double precision, parameter :: rhof_top = 1040.d0
  double precision, parameter :: rho_top = 2129.d0
  double precision, parameter :: rsm_top =25168.d0
  double precision, parameter :: alpha_top=0.58d0
  double precision, parameter :: rbM_top = 7.34d09
  double precision, parameter :: rlambdao_top =6.D08
  double precision, parameter :: rlambdac_top =rlambdao_top+alpha_top**2*rbM_top
  double precision, parameter :: ro11_t=rho_top+phi_top*rhof_top*(a_top-2.d0)
  double precision, parameter :: ro12_t=phi_top*rhof_top*(1.d0-a_top)
  double precision, parameter :: ro22_t=a_top*phi_top*rhof_top
  double precision, parameter :: lambda_t=rlambdao_top+rbM_top*(alpha_top-phi_top)**2
  double precision, parameter :: R_t=rbM_top*phi_top**2
  double precision, parameter :: ga_t=rbM_top*phi_top*(alpha_top-phi_top)
  double precision, parameter :: S_t=lambda_t+2*rmu_top
  double precision, parameter :: c1_t=S_t*R_t-ga_t**2
  double precision, parameter :: b1_t=-S_t*ro22_t-R_t*ro11_t+2*ga_t*ro12_t
  double precision, parameter :: a1_t=ro11_t*ro22_t-ro12_t**2
  double precision, parameter :: delta_t=b1_t**2-4.d0*a1_t*c1_t

  double precision:: cp_top
  double precision:: cps_top
  double precision:: cs_top

! total number of time steps
  integer, parameter :: NSTEP = 100000

! time step in seconds
  double precision, parameter :: DELTAT = 1.d-04

! parameters for the source
  double precision, parameter :: f0 = 80.d0
  double precision, parameter :: t0 = 1.d0/f0
  double precision, parameter :: factor =1.d02

! source
  integer, parameter :: ISOURCE = NX/2+1
  integer, parameter :: JSOURCE = NY/2 +1
  integer, parameter :: IDEB =  NX / 2 + 1
  integer, parameter :: JDEB =  NY / 2 + 1
  double precision, parameter :: xsource = DELTAX * ISOURCE
  double precision, parameter :: ysource = DELTAY * JSOURCE
! angle of source force clockwise with respect to vertical (Y) axis
  double precision, parameter :: ANGLE_FORCE = 0.d0

! receivers
  integer, parameter :: NREC = 2
  double precision, parameter :: ydeb = NPOINTS_PML*DELTAY+10.D0   ! first receiver x in meters
  double precision, parameter :: yfin = NY*DELTAY-NPOINTS_PML*DELTAY-10.d0   ! first receiver x in meters
  double precision, parameter :: xdeb =xsource  ! first receiver y in meters
  double precision, parameter :: xfin =xdeb   ! first receiver y in meters

! display information on the screen from time to time
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

! main arrays
  double precision, dimension(0:NX+1,0:NY+1) :: vx,vy,sigmaxx,sigma2,alp_sigma2,sigmayy,sigmaxy,vnorm
  double precision, dimension(0:NX+1,0:NY+1) :: vxf,vyf
  double precision, dimension(0:NX+1,0:NY+1) :: rho,rhof,rsm,rmu,rlambdac,rbM,alpha,etaokappa,rlambdao

! to interpolate material parameters at the right location in the staggered grid cell
  double precision rho_half_x_half_y,rhof_half_x_half_y,rsm_half_x_half_y,etaokappa_half_x_half_y

! for evolution of total energy in the medium
  double precision epsilon_xx,epsilon_yy,epsilon_xy
  double precision, dimension(NSTEP) :: total_energy_kinetic,total_energy_potential
  double precision c33_half_y

! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
  double precision, parameter :: K_MAX_PML = 1.d0
  double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0)   ! from Festa and Vilotte

! 2D arrays for the memory variables
  double precision, dimension(0:NX+1,0:NY+1) :: gamma11,gamma22
  double precision, dimension(0:NX+1,0:NY+1) :: gamma12_1
  double precision, dimension(0:NX+1,0:NY+1) :: xi_1,xi_2

  double precision, dimension(0:NX+1,0:NY+1) :: &
     memory_dx_vx1,memory_dx_vx2,memory_dy_vx,memory_dx_vy,memory_dy_vy1,memory_dy_vy2, &
     memory_dx_sigmaxx,memory_dx_sigmayy,memory_dx_sigmaxy, &
     memory_dx_sigma2vx,memory_dx_sigma2vxf,memory_dy_sigma2vy,memory_dy_sigma2vyf, &
     memory_dy_sigmaxx,memory_dy_sigmayy,memory_dy_sigmaxy

! 1D arrays for the damping profiles
  double precision, dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half_x,K_x_half_x,alpha_x_half_x,a_x_half_x,b_x_half_x
  double precision, dimension(NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half_y,K_y_half_y,alpha_y_half_y,a_y_half_y,b_y_half_y

  double precision thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
  double precision Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized
  double precision value_dx_vx1,value_dx_vx2,value_dx_vy,value_dx_sigmaxx,value_dx_sigmaxy
  double precision value_dy_vy1,value_dy_vy2,value_dy_vx,value_dy_sigmaxx,value_dy_sigmaxy
  double precision value_dx_sigma2vxf,value_dy_sigma2vyf

! for the source
  double precision a,t,source_term

! for receivers
  double precision xspacerec,yspacerec,distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec
  double precision, dimension(NREC) :: xrec,yrec

! for seismograms
  double precision, dimension(NSTEP,NREC) :: sisvx,sisvy,sisp

  integer i,j,it,irec

  double precision Courant_number_bottom,Courant_number_top,velocnorm_all,max_amplitude
  double precision Dispersion_number_bottom,Dispersion_number_top

!---
!--- program starts here
!---
  cp_bottom=(-b1_b+sqrt(delta_b))/(2.d0*a1_b);
  cps_bottom=(-b1_b-sqrt(delta_b))/(2.d0*a1_b);
  cp_bottom=sqrt(cp_bottom)
  cps_bottom=sqrt(cps_bottom)
  cs_bottom=sqrt(rmu_bottom/(ro11_b-ro12_b**2/ro22_b))

  cp_top=(-b1_t+sqrt(delta_t))/(2.d0*a1_t);
  cps_top=(-b1_t-sqrt(delta_t))/(2.d0*a1_t);
  cp_top=sqrt(cp_top)
  cps_top=sqrt(cps_top)
  cs_top=sqrt(rmu_top/(ro11_t-ro12_t**2/ro22_t))

  print *,'cp_bottom= ',cp_bottom
  print *,'cps_bottom=',cps_bottom
  print *,'cs_bottom= ',cs_bottom
  print *,'cp_top= ',cp_top
  print *,'cps_top=',cps_top
  print *,'cs_top= ',cs_top

  print *,'rho_bottom= ',rho_bottom
  print *,'rsm_bottom= ',rsm_bottom
  print *,'rho_top= ',rho_top
  print *,'rsm_top= ',rsm_top
  print *,'rmu_bottom= ',rmu_bottom
  print *,'rlambdac_bottom= ',rlambdac_bottom
  print *,'rlambdao_bottom= ',rlambdao_bottom
  print *,'alpha_bottom= ',alpha_bottom
  print *,'rbM_bottom= ',rbM_bottom
  print *,'etaokappa_bottom= ',etaokappa_bottom
  print *,'rmu_top= ',rmu_top
  print *,'rlambdac_top= ',rlambdac_top
  print *,'rlambdao_top= ',rlambdao_top
  print *,'alpha_top= ',alpha_top
  print *,'rbM_top= ',rbM_top
  print *,'etaokappa_top= ',etaokappa_top

  print *, 'DELTAT CPML=', DELTAT
  print *,'2D poroelastic finite-difference code in velocity and stress formulation with C-PML'
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
  Rcoef = 0.001d0

! check that NPOWER is okay
  if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  if (HETEROGENEOUS_MODEL) then
    d0_x = - (NPOWER + 1) * max(cp_bottom,cp_top) * log(Rcoef) / (2.d0 * thickness_PML_x)
    d0_y = - (NPOWER + 1) * max(cp_bottom,cp_top) * log(Rcoef) / (2.d0 * thickness_PML_y)
  else
    d0_x = - (NPOWER + 1) * cp_bottom * log(Rcoef) / (2.d0 * thickness_PML_x)
    d0_y = - (NPOWER + 1) * cp_bottom * log(Rcoef) / (2.d0 * thickness_PML_y)
  endif

  print *,'d0_x = ',d0_x
  print *,'d0_y = ',d0_y

  d_x(:) = ZERO
  d_x_half_x(:) = ZERO

  d_y(:) = ZERO
  d_y_half_y(:) = ZERO

  K_x(:) = 1.d0
  K_x_half_x(:) = 1.d0

  K_y(:) = 1.d0
  K_y_half_y(:) = 1.d0

  alpha_x(:) = ZERO
  alpha_x_half_x(:) = ZERO

  alpha_y(:) = ZERO
  alpha_y_half_y(:) = ZERO

  a_x(:) = ZERO
  a_x_half_x(:) = ZERO

  a_y(:) = ZERO
  a_y_half_y(:) = ZERO

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x
  xoriginright = (NX-1)*DELTAX - thickness_PML_x

  do i = 1,NX

! abscissa of current grid point along the damping profile
    xval = DELTAX * dble(i-1)

!!!! ---------- left edge
    if (USE_PML_LEFT) then

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
        d_x_half_x(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!!!! ---------- right edge
    if (USE_PML_RIGHT) then

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
        d_x_half_x(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_x_half_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

! just in case, for -5 at the end
    if (alpha_x(i) < ZERO) alpha_x(i) = ZERO
    if (alpha_x_half_x(i) < ZERO) alpha_x_half_x(i) = ZERO

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
    b_x_half_x(i) = exp(- (d_x_half_x(i) / K_x_half_x(i) + alpha_x_half_x(i)) * DELTAT)

! this to avoid division by zero outside the PML
    if (abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) /&
      (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if (abs(d_x_half_x(i)) > 1.d-6) a_x_half_x(i) = d_x_half_x(i)&
     * (b_x_half_x(i) - 1.d0) / (K_x_half_x(i) * (d_x_half_x(i) + K_x_half_x(i) * alpha_x_half_x(i)))

  enddo

!!!!!!!!!!!!! added Y damping profile

! origin of the PML layer (position of right edge minus thickness, in meters)
  yoriginbottom = thickness_PML_y
  yorigintop = NY*DELTAY - thickness_PML_y

  do j = 1,NY

! abscissa of current grid point along the damping profile
    yval = DELTAY * dble(j-1)

!!!! ---------- bottom edge
    if (USE_PML_BOTTOM) then

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
        d_y_half_y(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!!!! ---------- top edge
    if (USE_PML_TOP) then

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
        d_y_half_y(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
        K_y_half_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

! just in case, for -5 at the end
!   if (alpha_y(j) < ZERO) alpha_y(j) = ZERO
!   if (alpha_y_half_y(j) < ZERO) alpha_y_half_y(j) = ZERO

    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT)
    b_y_half_y(j) = exp(- (d_y_half_y(j) / K_y_half_y(j) + alpha_y_half_y(j)) * DELTAT)

! this to avoid division by zero outside the PML
    if (abs(d_y(j)) > 1.d-6) a_y(j) = d_y(j) * (b_y(j) - 1.d0) &
     / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
    if (abs(d_y_half_y(j)) > 1.d-6) a_y_half_y(j) = d_y_half_y(j)&
      * (b_y_half_y(j) - 1.d0) / (K_y_half_y(j) * (d_y_half_y(j) + K_y_half_y(j) * alpha_y_half_y(j)))

  enddo

! compute the Lame parameters and density
  do j = 0,NY+1
    do i = 0,NX+1
      if (HETEROGENEOUS_MODEL .and. DELTAY*dble(j-1) > INTERFACE_HEIGHT) then
         rho(i,j)= rho_top
         rhof(i,j) = rhof_top
         rsm(i,j) = rsm_top
         rmu(i,j)= rmu_top
         rlambdac(i,j) = rlambdac_top
         rbM(i,j) = rbM_top
         alpha(i,j)=alpha_top
         etaokappa(i,j)=etaokappa_top
         rlambdao(i,j) = rlambdao_top
      else
         rho(i,j)= rho_bottom
         rhof(i,j) = rhof_bottom
         rsm(i,j) = rsm_bottom
         rmu(i,j)= rmu_bottom
         rlambdac(i,j) = rlambdac_bottom
         rbM(i,j) = rbM_bottom
         alpha(i,j)=alpha_bottom
         etaokappa(i,j)=etaokappa_bottom
         rlambdao(i,j) = rlambdao_bottom
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
  enddo

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
  Courant_number_bottom = cp_bottom * DELTAT / min(DELTAX,DELTAY)
  Dispersion_number_bottom=min(cs_bottom,cps_bottom)/(2.5d0*f0*max(DELTAX,DELTAY))
  print *,'Courant number at the bottom is ',Courant_number_bottom
  print *,'Dispersion number at the bottom is ',Dispersion_number_bottom
  print *
  if (Courant_number_bottom > 1.d0/sqrt(2.d0)) stop 'time step is too large, simulation will be unstable'

  if (HETEROGENEOUS_MODEL) then
    Courant_number_top = max(cp_top,cp_bottom) * DELTAT / min(DELTAX,DELTAY)
    Dispersion_number_top=min(cs_top,cs_bottom,cps_bottom,cps_top)/(2.5d0*f0*max(DELTAX,DELTAY))
    print *,'Courant number at the top is ',Courant_number_top
    print *
    print *,'Dispersion number at the top is ',Dispersion_number_top
    if (Courant_number_top > 6.d0/7.d0/sqrt(2.d0)) stop 'time step is too large, simulation will be unstable'
  endif

! suppress old files
! call system('rm -f Vx_*.dat Vy_*.dat vect*.ps image*.pnm image*.gif')

! initialize arrays
  vx(:,:) = ZERO
  vy(:,:) = ZERO
  sigmaxx(:,:) = ZERO
  sigmayy(:,:) = ZERO
  sigmaxy(:,:) = ZERO
  sigma2(:,:) = ZERO
  alp_sigma2(:,:) = ZERO
  gamma11(:,:)=0.d0
  gamma22(:,:)=0.d0
  gamma12_1(:,:)=0.d0
  gamma12_1(:,:)=0.d0
  xi_1(:,:)=0.d0
  xi_2(:,:)=0.d0
  vxf(:,:) = ZERO
  vyf(:,:) = ZERO

     memory_dx_vx1(:,:)=0.d0
     memory_dx_vx2(:,:)=0.d0
     memory_dy_vx(:,:)=0.d0
     memory_dx_vy(:,:)=0.d0
     memory_dy_vy1(:,:)=0.d0
     memory_dy_vy2(:,:)=0.d0
     memory_dx_sigmaxx(:,:)=0.d0
     memory_dx_sigmayy(:,:)=0.d0
     memory_dx_sigmaxy(:,:)=0.d0
     memory_dx_sigma2vx(:,:)=0.d0
     memory_dx_sigma2vxf(:,:)=0.d0
     memory_dy_sigmaxx(:,:)=0.d0
     memory_dy_sigmayy(:,:)=0.d0
     memory_dy_sigmaxy(:,:)=0.d0
     memory_dy_sigma2vy(:,:)=0.d0
     memory_dy_sigma2vyf(:,:)=0.d0

! initialize seismograms
  sisvx(:,:) = ZERO
  sisvy(:,:) = ZERO
  sisp(:,:) = ZERO

! initialize total energy
  total_energy_kinetic(:) = ZERO
  total_energy_potential(:) = ZERO

!---
!---  beginning of time loop
!---

  do it = 1,NSTEP

!----------------------
! compute stress sigma
!----------------------

!-----------------------------------
! update memory variables for C-PML
!-----------------------------------

  do j = 2,NY
    do i = 1,NX-1

!  memory of sigmaxx
      value_dx_sigmaxx =(27.d0*vx(i+1,j)-27.d0*vx(i,j)-vx(i+2,j)+vx(i-1,j))/DELTAX/24.D0
      value_dy_sigmaxx =(27.d0*vy(i,j)-27.d0*vy(i,j-1)-vy(i,j+1)+vy(i,j-2))/DELTAY/24.D0

      memory_dx_sigmaxx(i,j) = b_x_half_x(i) * memory_dx_sigmaxx(i,j) + a_x_half_x(i) * value_dx_sigmaxx
      memory_dy_sigmaxx(i,j) = b_y(j) * memory_dy_sigmaxx(i,j) + a_y(j) * value_dy_sigmaxx


      gamma11(i,j) = gamma11(i,j)+DELTAT*(value_dx_sigmaxx / K_x_half_x(i) + memory_dx_sigmaxx(i,j))

      gamma22(i,j) = gamma22(i,j)+DELTAT*(value_dy_sigmaxx / K_y(j) + memory_dy_sigmaxx(i,j))

! sigma2
      value_dx_sigma2vxf=(27.d0*vxf(i+1,j)-27.d0* vxf(i,j)-vxf(i+2,j)+vxf(i-1,j)) / DELTAX/24.D0
      value_dy_sigma2vyf=(27.d0*vyf(i,j)-27.d0*vyf(i,j-1)-vyf(i,j+1)+vyf(i,j-2)) / DELTAY/24.D0

      memory_dx_sigma2vxf(i,j) = b_x_half_x(i) * memory_dx_sigma2vxf(i,j) + a_x_half_x(i) * value_dx_sigma2vxf
      memory_dy_sigma2vyf(i,j) = b_y(j) * memory_dy_sigma2vyf(i,j) + a_y(j) * value_dy_sigma2vyf

      xi_1(i,j) = xi_1(i,j) -(value_dx_sigma2vxf/ K_x_half_x(i) + memory_dx_sigma2vxf(i,j))*DELTAT

      xi_2(i,j) = xi_2(i,j) -(value_dy_sigma2vyf/K_y(j)+memory_dy_sigma2vyf(i,j))*DELTAT

    sigma2(i,j)=-alpha(i,j)*rbM(i,j)*(gamma11(i,j)+gamma22(i,j))+rbM(i,j)*(xi_1(i,j)+xi_2(i,j))

    enddo
  enddo

! add the source (point source located at a given grid point)
  a = pi*pi*f0*f0
  t = dble(it-1)*DELTAT

! Gaussian
  source_term = factor * exp(-a*(t-t0)**2)/(-2.d0*a)

! first derivative of a Gaussian
! source_term =  factor * 2.d0*a*(t-t0)*exp(-a*(t-t0)**2)
! source_term =  factor *(t-t0)*exp(-a*(t-t0)**2)

! Ricker source time function (second derivative of a Gaussian)
! source_term = factor * (1.d0 - 2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)

! define location of the source
  i = ISOURCE
  j = JSOURCE

! add the source term
  sigma2(i,j) = sigma2(i,j) + source_term*rbM(i,j)

  do j = 1,NY-1
    do i = 2,NX

! interpolate material parameters at the right location in the staggered grid cell
      c33_half_y = 2.d0/(1.d0/rmu(i,j)+1.d0/rmu(i,j+1))
      c33_half_y = rmu(i,j+1)

      value_dx_sigmaxy = (27.d0*vy(i,j) - 27.d0*vy(i-1,j)-vy(i+1,j)+vy(i-2,j)) / DELTAX/24.D0
      value_dy_sigmaxy = (27.d0*vx(i,j+1) - 27.d0*vx(i,j)-vx(i,j+2)+vx(i,j-1)) / DELTAY/24.D0

      memory_dx_sigmaxy(i,j) = b_x(i) * memory_dx_sigmaxy(i,j) + a_x(i) * value_dx_sigmaxy
      memory_dy_sigmaxy(i,j) = b_y_half_y(j) * memory_dy_sigmaxy(i,j) + a_y_half_y(j) * value_dy_sigmaxy

      sigmaxy(i,j) = sigmaxy(i,j) + &
      c33_half_y/1.d0 * (value_dx_sigmaxy / K_x(i) + memory_dx_sigmaxy(i,j) + &
        value_dy_sigmaxy / K_y(j) + memory_dy_sigmaxy(i,j)) * DELTAT

    enddo
  enddo

 do j = 2,NY
    do i = 1,NX-1
      sigmaxx(i,j)=(rlambdao(i,j)+2.d0*rmu(i,j))*gamma11(i,j)+rlambdao(i,j)*gamma22(i,j) -alpha(i,j)*sigma2(i,j)
      sigmayy(i,j)=rlambdao(i,j)*gamma11(i,j)+(rlambdao(i,j)+2.d0*rmu(i,j))*gamma22(i,j) -alpha(i,j)*sigma2(i,j)
    enddo
  enddo

!------------------
! compute velocity
!------------------

!-----------------------------------
! update memory variables for C-PML
!-----------------------------------

  do j = 2,NY
    do i = 2,NX
    co=(rho(i,j)*rsm(i,j)-rhof(i,j)*rhof(i,j))/DELTAT
    c1=co+rho(i,j)*etaokappa(i,j)*0.5d0
    c2=co-rho(i,j)*etaokappa(i,j)*0.5d0
    vtemp=vxf(i,j)
      value_dx_vx1 = (27.d0*sigmaxx(i,j) - 27.d0*sigmaxx(i-1,j)&
      -sigmaxx(i+1,j)+sigmaxx(i-2,j)) / DELTAX/24.D0
      value_dx_vx2 = (27.d0*sigma2(i,j) - 27.d0*sigma2(i-1,j)-sigma2(i+1,j)+sigma2(i-2,j)) / DELTAX/24.D0
      value_dy_vx = (27.d0*sigmaxy(i,j) - 27.d0*sigmaxy(i,j-1)-sigmaxy(i,j+1)+sigmaxy(i,j-2)) / DELTAY/24.D0

      memory_dx_vx1(i,j) = b_x(i) * memory_dx_vx1(i,j) + a_x(i) * value_dx_vx1
      memory_dx_vx2(i,j) = b_x(i) * memory_dx_vx2(i,j) + a_x(i) * value_dx_vx2
      memory_dy_vx(i,j) = b_y(j) * memory_dy_vx(i,j) + a_y(j) * value_dy_vx

      vxf(i,j) = (c2*vxf(i,j) + &
         (-rhof(i,j)*(value_dx_vx1/ K_x(i) + memory_dx_vx1(i,j) &
          + value_dy_vx / K_y(j) + memory_dy_vx(i,j)) &
          -rho(i,j)*(value_dx_vx2/ K_x(i) + memory_dx_vx2(i,j)) &
         )) /c1

      vtemp=(vtemp+vxf(i,j))*0.5d0

      vx(i,j) = vx(i,j) + &
         (rsm(i,j)*(value_dx_vx1/ K_x(i) + memory_dx_vx1(i,j)+ &
          value_dy_vx / K_y(j) + memory_dy_vx(i,j))+&
          rhof(i,j)*(value_dx_vx2/ K_x(i) + memory_dx_vx2(i,j)) + &
          rhof(i,j)*etaokappa(i,j)*vtemp)&
         /co

    enddo
  enddo

  do j = 1,NY-1
    do i = 1,NX-1

      rho_half_x_half_y = rho(i,j+1)
      rsm_half_x_half_y = rsm(i,j+1)
      rhof_half_x_half_y = rhof(i,j+1)
      etaokappa_half_x_half_y = etaokappa(i,j+1)

      co=(rho_half_x_half_y*rsm_half_x_half_y-rhof_half_x_half_y**2)/DELTAT
      c1=co+rho_half_x_half_y*etaokappa_half_x_half_y*0.5d0
      c2=co-rho_half_x_half_y*etaokappa_half_x_half_y*0.5d0
      vtemp=vyf(i,j)

      value_dx_vy = (27.d0*sigmaxy(i+1,j) - 27.d0*sigmaxy(i,j)-sigmaxy(i+2,j)+sigmaxy(i-1,j)) / DELTAX/24.D0
      value_dy_vy1 = (27.d0*sigmayy(i,j+1)- 27.d0*sigmayy(i,j)&
      -sigmayy(i,j+2)+sigmayy(i,j-1)) / DELTAY/24.D0
      value_dy_vy2 = (27.d0*sigma2(i,j+1) - 27.d0*sigma2(i,j)-sigma2(i,j+2)+sigma2(i,j-1)) / DELTAY/24.D0

      memory_dx_vy(i,j)  = b_x_half_x(i) * memory_dx_vy(i,j) + a_x_half_x(i) * value_dx_vy
      memory_dy_vy1(i,j) = b_y_half_y(j) * memory_dy_vy1(i,j) + a_y_half_y(j) * value_dy_vy1
      memory_dy_vy2(i,j) = b_y_half_y(j) * memory_dy_vy2(i,j) + a_y_half_y(j) * value_dy_vy2

   vyf(i,j) = (c2*vyf(i,j) + &
 (-rhof_half_x_half_y*(value_dx_vy / K_x_half_x(i) + memory_dx_vy(i,j) &
 +value_dy_vy1 / K_y_half_y(j) + memory_dy_vy1(i,j))&
  -rho_half_x_half_y*(value_dy_vy2 / K_y_half_y(j) + memory_dy_vy2(i,j)))&
  ) /c1
      vtemp=(vtemp+vyf(i,j))*0.5d0

   vy(i,j) = vy(i,j) + &
 (rsm_half_x_half_y*(value_dx_vy / K_x_half_x(i) + memory_dx_vy(i,j)&
+ value_dy_vy1 / K_y_half_y(j) + memory_dy_vy1(i,j))&
+ rhof_half_x_half_y*(value_dy_vy2 / K_y_half_y(j) + memory_dy_vy2(i,j))&
+ rhof_half_x_half_y*etaokappa_half_x_half_y*vtemp)&
 /co

    enddo
  enddo

! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
  vx(1,:) = ZERO
  vx(NX,:) = ZERO

  vx(:,1) = ZERO
  vx(:,NY) = ZERO

  vy(1,:) = ZERO
  vy(NX,:) = ZERO

  vy(:,1) = ZERO
  vy(:,NY) = ZERO

  vxf(1,:) = ZERO
  vxf(NX,:) = ZERO

  vxf(:,1) = ZERO
  vxf(:,NY) = ZERO

  vyf(1,:) = ZERO
  vyf(NX,:) = ZERO

  vyf(:,1) = ZERO
  vyf(:,NY) = ZERO

! store seismograms
  do irec = 1,NREC
! x component
    sisvx(it,irec) = vx(ix_rec(irec),iy_rec(irec))
! y component
    sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec))
! fluid pressure
    sisp(it,irec) = sigma2(ix_rec(irec),iy_rec(irec))
  enddo

! compute total energy

! compute kinetic energy first, defined as 1/2 rho ||v||^2
! in principle we should use rho_half_x_half_y instead of rho for vy
! in order to interpolate density at the right location in the staggered grid cell
! but in a homogeneous medium we can safely ignore it
 total_energy_kinetic(it) = 0.5d0 * &
sum(rho(NPOINTS_PML:NX-NPOINTS_PML+1,NPOINTS_PML:NY-NPOINTS_PML+1)&
*(vx(NPOINTS_PML:NX-NPOINTS_PML+1,NPOINTS_PML:NY-NPOINTS_PML+1)**2&
+vy(NPOINTS_PML:NX-NPOINTS_PML+1,NPOINTS_PML:NY-NPOINTS_PML+1)**2))&
*DELTAX * DELTAY+&
0.5d0*sum(rsm(NPOINTS_PML:NX-NPOINTS_PML+1,NPOINTS_PML:NY-NPOINTS_PML+1)&
*(vxf(NPOINTS_PML:NX-NPOINTS_PML+1,NPOINTS_PML:NY-NPOINTS_PML+1)**2&
+vyf(NPOINTS_PML:NX-NPOINTS_PML+1,NPOINTS_PML:NY-NPOINTS_PML+1)**2))&
*DELTAX*DELTAY

! add potential energy, defined as 1/2 epsilon_ij sigma_ij
! in principle we should interpolate the medium parameters at the right location
! in the staggered grid cell but in a homogeneous medium we can safely ignore it
  total_energy_potential(it) = ZERO

  do j = NPOINTS_PML,NY-NPOINTS_PML+1
    do i = NPOINTS_PML,NX-NPOINTS_PML+1
      epsilon_xx = ((rlambdao(i,j) + 2.d0*rmu(i,j)) * sigmaxx(i,j) - rlambdao(i,j) * sigmayy(i,j)) / &
        (4.d0 * rmu(i,j) * (rlambdao(i,j) + rmu(i,j)))
      epsilon_yy = ((rlambdao(i,j) + 2.d0*rmu(i,j)) * sigmayy(i,j) - rlambdao(i,j) * sigmaxx(i,j)) / &
        (4.d0 * rmu(i,j) * (rlambdao(i,j) + rmu(i,j)))
      epsilon_xy = sigmaxy(i,j) / (2.d0 * rmu(i,j))
      total_energy_potential(it) = total_energy_potential(it) + &
        0.5d0 * (epsilon_xx * sigmaxx(i,j) + epsilon_yy * sigmayy(i,j) + 2.d0 * epsilon_xy * sigmaxy(i,j)&
        +sigma2(i,j)**2/rbM(i,j)&
        +2.d0*rhof(i,j)*(vx(i,j)*vxf(i,j)+vy(i,j)*vyf(i,j)))*DELTAX * DELTAY
    enddo
  enddo

! output information
  if (mod(it,IT_DISPLAY) == 0 .or. it == 5) then

! print maximum of norm of velocity
    velocnorm_all = maxval(sqrt(vx(:,:)**2 + vy(:,:)**2))
    print *,'Time step # ',it,' out of ',NSTEP
    print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
    print *,'Max norm velocity vector V (m/s) = ',velocnorm_all
    print *,'total energy = ',total_energy_kinetic(it) + total_energy_potential(it)
    print *

! check stability of the code, exit if unstable
    if (velocnorm_all > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

    vnorm(:,:)=sqrt(vx(:,:)**2+vy(:,:)**2)

  call create_color_image(vx,NX+2,NY+2,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
  NPOINTS_PML,USE_PML_LEFT,USE_PML_RIGHT,USE_PML_BOTTOM, &
  USE_PML_TOP,1,max_amplitude,JINTERFACE)

  call create_color_image(vy,NX+2,NY+2,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
  NPOINTS_PML,USE_PML_LEFT,USE_PML_RIGHT,USE_PML_BOTTOM, &
  USE_PML_TOP,2,max_amplitude,JINTERFACE)

! save temporary partial seismograms to monitor the behavior of the simulation
! while it is running
  call write_seismograms(sisvx,sisvy,sisp,NSTEP,NREC,DELTAT,t0)

  endif

  enddo   ! end of time loop

! save seismograms
  call write_seismograms(sisvx,sisvy,sisp,NSTEP,NREC,DELTAT,t0)

! save total energy
  open(unit=20,file='energy.dat',status='unknown')
  do it = 1,NSTEP
    write(20,*) sngl(dble(it-1)*DELTAT), sngl(total_energy_kinetic(it) + total_energy_potential(it))
  enddo
  close(20)

! create script for Gnuplot for total energy
  open(unit=20,file='plot_energy',status='unknown')
  write(20,*) 'set term x11'
  write(20,*) '# set term postscript landscape monochrome dashed "Helvetica" 22'
  write(20,*) '# set xrange [0:7]'
  write(20,*) '# set yrange [-4:4.5]'
  write(20,*)
  write(20,*) 'set xlabel "Time (s)"'
  write(20,*) 'set ylabel "Total energy"'
  write(20,*)
  write(20,*) '# set output "cpml_total_energy.eps"'
  write(20,*) 'plot "energy.dat" us 1:2 t ''Ec'' w l lc 1, "energy.dat" us 1:3 &
    & t ''Ep'' w l lc 3, "energy.dat" us 1:4 t ''Total energy'' w l lc 4'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)
  close(20)

! create script for Gnuplot
  open(unit=20,file='plotgnu',status='unknown')
  write(20,*) 'set term x11'
  write(20,*) '# set term postscript landscape monochrome dashed "Helvetica" 22'
  write(20,*) '#set xrange [0:7]'
  write(20,*)
  write(20,*) 'set xlabel "Time (s)"'
  write(20,*) 'set ylabel "Amplitude (m / s)"'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vx_receiver_001.eps"'
  write(20,*) '#set yrange [-4:4.5]'
  write(20,*) 'plot "Vx_file_001.dat" t ''Vx C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_001.eps"'
  write(20,*) '#set yrange [-15:19]'
  write(20,*) 'plot "Vy_file_001.dat" t ''Vy C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vx_receiver_002.eps"'
  write(20,*) '#set yrange [-12:16]'
  write(20,*) 'plot "Vx_file_002.dat" t ''Vx C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  write(20,*) 'set output "v_sigma_Vy_receiver_002.eps"'
  write(20,*) '#set yrange [-7:10]'
  write(20,*) 'plot "Vy_file_002.dat" t ''Vy C-PML'' w l lc 1'
  write(20,*) 'pause -1 "Hit any key..."'
  write(20,*)

  close(20)

  print *
  print *,'End of the simulation'
  print *

  end program seismic_CPML_2D_poroelastic_fourth


!----
!----  save the seismograms in ASCII text format
!----

  subroutine write_seismograms(sisvx,sisvy,sisp,nt,nrec,DELTAT,t0)

  implicit none

  integer nt,nrec
  double precision DELTAT,t0

  double precision sisvx(nt,nrec)
  double precision sisvy(nt,nrec)
  double precision sisp(nt,nrec)

  integer irec,it

  character(len=100) file_name

! X component
  do irec=1,nrec
    write(file_name,"('Vx_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvx(it,irec))
    enddo
    close(11)
  enddo

! Z component
  do irec=1,nrec
    write(file_name,"('Vy_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvy(it,irec))
    enddo
    close(11)
  enddo

! fluid pressure
  do irec=1,nrec
    write(file_name,"('Pf_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisp(it,irec))
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
  double precision, parameter :: POWER_DISPLAY = 0.25d0

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
    write(file_name,"('image',i5.5,'_Vx.pnm')") it
    write(system_command,"('convert image',i5.5,'_Vx.pnm image',i5.5,'_Vx.gif ; rm image',i5.5,'_Vx.pnm')") it,it,it
  endif
  if (field_number == 2) then
    write(file_name,"('image',i5.5,'_Vy.pnm')") it
    write(system_command,"('convert image',i5.5,'_Vy.pnm image',i5.5,'_Vy.gif ; rm image',i5.5,'_Vy.pnm')") it,it,it
  endif
  if (field_number == 3) then
    write(file_name,"('image',i5.5,'_Vnorm.pnm')") it
    write(system_command,"('convert image',i5.5,'_Vnorm.pnm image',i5.5,'_Vnorm.gif ; rm image',i5.5,'_Vnorm.pnm')") it,it,it
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

