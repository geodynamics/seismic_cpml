!
! SEISMIC_CPML Version 1.1.3, September 2015.
!
! Copyright CNRS, France.
! Contributor: Dimitri Komatitsch, komatitsch aT lma DOT cnrs-mrs DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional heterogeneous acoustic wave equation
! using a finite-difference method with Convolutional Perfectly Matched
! Layer (C-PML) conditions.
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  program seismic_CPML_2D_pressure

! 2D acoustic finite-difference code in pressure formulation
! with Convolutional-PML (C-PML) absorbing conditions for an heterogeneous acoustic medium

! Dimitri Komatitsch, CNRS, Marseille, September 2015.

! The pressure wave equation in an inviscid heterogeneous fluid is:
!
! 1/Kappa d2p / dt2 = div(grad(p) / rho) = d(1/rho dp/dx)/dx + d(1/rho dp/dy)/dy
!
! (see for instance Komatitsch and Tromp, Geophysical Journal International, vol. 149, p. 390-412 (2002), equations (19) and (21))
!
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
!            |                   |
!      dp/dy +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!            p       dp/dx

! pour afficher le resultat 2D sous forme d'image en couleurs :
!
!   " display image*.gif " ou " gimp image*.gif "
!
!   " montage -geometry +0+3 -tile 1x60 image*.gif allfiles.gif "
!   puis " display allfiles.gif " ou bien " gimp allfiles.gif "
!
! When compiling with Intel ifort, use " -assume byterecl " option to create temporary binary PNM images

  implicit none

! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.

! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 10

! total number of grid points in each direction of the grid
  integer, parameter :: NX = 440
  integer, parameter :: NY = 440

! size of a grid cell
  double precision, parameter :: DELTAX = 12.5d0, DELTAY = DELTAX

! P-velocity, S-velocity and density
  double precision, parameter :: cp_value = 2500.d0
  double precision, parameter :: rho_value = 2200.d0

! total number of time steps
  integer, parameter :: NSTEP = 1000

! time step in seconds
  double precision, parameter :: DELTAT = 2.d-3

! rapport DELTAT / DELTAX
  double precision, parameter :: DELTAT_OVER_DELTAX = DELTAT / DELTAX

! parameters for the source
! tres peu de dispersion si 8 Hz, correct si 10 Hz, dispersion sur S si 12 Hz
! mauvais si 16 Hz, tres mauvais si 25 Hz
  double precision, parameter :: f0 = 8.d0
! decaler le Ricker en temps pour s'assurer qu'il est nul a t = 0
  double precision, parameter :: t0 = 1.2d0 / f0
  double precision, parameter :: factor = 1.d4

! source in the medium
  double precision, parameter :: XSOURCE = (NX / 2) * DELTAX
  double precision, parameter :: YSOURCE = (NY / 2) * DELTAX

! receivers
  integer, parameter :: NREC = 11
  double precision, parameter :: xdeb = 100.d0   ! first receiver x in meters
  double precision, parameter :: ydeb = 2400.d0  ! first receiver y in meters
  double precision, parameter :: xfin = 3850.d0  ! last receiver x in meters
  double precision, parameter :: yfin = 2400.d0  ! last receiver y in meters

! affichage periodique d'informations a l'ecran
  integer, parameter :: IT_AFFICHE = 100

! zero
  double precision, parameter :: ZERO = 0.d0

! large value for maximum
  double precision, parameter :: HUGEVAL = 1.d+30

! pressure threshold above which we consider the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! valeur de PI
  double precision, parameter :: PI = 3.141592653589793238462643d0

! main arrays
  double precision, dimension(NX,NY) :: xgrid,ygrid, &
      pressurepast,pressurepresent,pressurefuture, &
      pressure_xx,pressure_yy,dpressurexx_dx,dpressureyy_dy,kappa,cp,rho

! for seismograms
  double precision, dimension(NSTEP,NREC) :: sispressure

! for source
  integer i_source,j_source
  double precision a,t,force_x,force_y

! for receivers
  double precision xspacerec,yspacerec,distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec
  double precision, dimension(NREC) :: xrec,yrec

  integer i,j,it,irec

  double precision nombre_Courant,pressure_max_all

! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

  double precision, parameter :: K_MAX_PML = 1.d0 ! from Gedney page 8.11
  double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from festa and Vilotte

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
  double precision, dimension(NX,NY) :: &
      memory_dpressure_dx, &
      memory_dpressure_dy, &
      memory_dpressurexx_dx, &
      memory_dpressureyy_dy

  double precision :: &
      value_dpressure_dx, &
      value_dpressure_dy, &
      value_dpressurexx_dx, &
      value_dpressureyy_dy

! 1D arrays for the damping profiles
  double precision, dimension(NX) :: d_x,K_x,alpha_prime_x,a_x,b_x,d_x_half,K_x_half,alpha_prime_x_half,a_x_half,b_x_half
  double precision, dimension(NY) :: d_y,K_y,alpha_prime_y,a_y,b_y,d_y_half,K_y_half,alpha_prime_y_half,a_y_half,b_y_half

  double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
  double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

!---
!--- program starts here
!---

  print *
  print *,'2D acoustic finite-difference code in pressure formulation'
  print *

! create the mesh of grid points
  do j = 1,NY
    do i = 1,NX
      xgrid(i,j) = DELTAX * dble(i-1)
      ygrid(i,j) = DELTAY * dble(j-1)
    enddo
  enddo

! compute the Lame parameters and density
  do j = 1,NY
    do i = 1,NX

! one can change the values of the density and P velocity model here to make it heterogeneous
      cp(i,j) = cp_value
      rho(i,j) = rho_value

      kappa(i,j) = rho(i,j)*cp(i,j)*cp(i,j)
    enddo
  enddo

! find closest grid point for the source
  dist = HUGEVAL
  do j=1,NY
    do i=1,NX
      distval = sqrt((xgrid(i,j)-XSOURCE)**2 + (ygrid(i,j)-YSOURCE)**2)
      if(distval < dist) then
        dist = distval
        i_source = i
        j_source = j
      endif
    enddo
  enddo
  print *,'closest grid point for the source found at distance ',dist,' in i,j = ',i_source,j_source

! define location of receivers
  print *
  print *,'There are ',NREC,' receivers'
  print *
  xspacerec = (xfin-xdeb) / dble(NREC-1)
  yspacerec = (yfin-ydeb) / dble(NREC-1)
  do irec=1,NREC
    xrec(irec) = xdeb + dble(irec-1)*xspacerec
    yrec(irec) = ydeb + dble(irec-1)*yspacerec
  enddo

! find closest grid point for each receiver
  do irec=1,NREC
    dist = HUGEVAL
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((xgrid(i,j)-xrec(irec))**2 + (ygrid(i,j)-yrec(irec))**2)
      if(distval < dist) then
        dist = distval
        ix_rec(irec) = i
        iy_rec(irec) = j
      endif
    enddo
    enddo
    print *,'closest grid point for receiver ',irec,' found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec)
  enddo

! afficher la taille du modele
  print *
  print *,'taille du modele suivant X = ',maxval(xgrid)
  print *,'taille du modele suivant Y = ',maxval(ygrid)
  print *

! verifier que la condition de stabilite de Courant est respectee
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
  nombre_Courant = maxval(cp) * DELTAT_OVER_DELTAX
  print *,'le nombre de Courant vaut ',nombre_Courant
  print *
  if(nombre_Courant > 1.d0/sqrt(2.d0)) then
    stop 'le pas de temps est trop grand, simulation instable'
  endif

! initialiser les tableaux
  pressurepresent(:,:) = ZERO
  pressurepast(:,:) = ZERO

! initialiser sismogrammes
  sispressure(:,:) = ZERO

! PML
  memory_dpressure_dx(:,:) = ZERO
  memory_dpressure_dy(:,:) = ZERO
  memory_dpressurexx_dx(:,:) = ZERO
  memory_dpressureyy_dy(:,:) = ZERO

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
  thickness_PML_x = NPOINTS_PML * DELTAX
  thickness_PML_y = NPOINTS_PML * DELTAY

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0

! check that NPOWER is okay
  if(NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * maxval(cp) * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * maxval(cp) * log(Rcoef) / (2.d0 * thickness_PML_y)

  print *,'d0_x = ',d0_x
  print *,'d0_y = ',d0_y
  print *

  d_x(:) = ZERO
  d_x_half(:) = ZERO
  K_x(:) = 1.d0
  K_x_half(:) = 1.d0
  alpha_prime_x(:) = ZERO
  alpha_prime_x_half(:) = ZERO
  a_x(:) = ZERO
  a_x_half(:) = ZERO

  d_y(:) = ZERO
  d_y_half(:) = ZERO
  K_y(:) = 1.d0
  K_y_half(:) = 1.d0
  alpha_prime_y(:) = ZERO
  alpha_prime_y_half(:) = ZERO
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
    if(USE_PML_XMIN) then

! define damping profile at the grid points
      abscissa_in_PML = xoriginleft - xval
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0)
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- right edge
    if(USE_PML_XMAX) then

! define damping profile at the grid points
      abscissa_in_PML = xval - xoriginright
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

! just in case, for -5 at the end
    if(alpha_prime_x(i) < ZERO) alpha_prime_x(i) = ZERO
    if(alpha_prime_x_half(i) < ZERO) alpha_prime_x_half(i) = ZERO

    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_prime_x(i)) * DELTAT)
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_prime_x_half(i)) * DELTAT)

! this to avoid division by zero outside the PML
    if(abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_prime_x(i)))
    if(abs(d_x_half(i)) > 1.d-6) a_x_half(i) = d_x_half(i) * &
      (b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_prime_x_half(i)))

  enddo

! damping in the Y direction

! origin of the PML layer (position of right edge minus thickness, in meters)
  yoriginbottom = thickness_PML_y
  yorigintop = (NY-1)*DELTAY - thickness_PML_y

  do j = 1,NY

! abscissa of current grid point along the damping profile
    yval = DELTAY * dble(j-1)

!---------- bottom edge
    if(USE_PML_YMIN) then

! define damping profile at the grid points
      abscissa_in_PML = yoriginbottom - yval
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0)
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- top edge
    if(USE_PML_YMAX) then

! define damping profile at the grid points
      abscissa_in_PML = yval - yorigintop
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! define damping profile at half the grid points
      abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop
      if(abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! this taken from Gedney page 8.2
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_prime_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_prime_y(j)) * DELTAT)
    b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_prime_y_half(j)) * DELTAT)

! this to avoid division by zero outside the PML
    if(abs(d_y(j)) > 1.d-6) a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_prime_y(j)))
    if(abs(d_y_half(j)) > 1.d-6) a_y_half(j) = d_y_half(j) * &
      (b_y_half(j) - 1.d0) / (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_prime_y_half(j)))

  enddo

! beginning of the time loop
  do it = 1,NSTEP

! calculer les premieres derivees spatiales multipliees par les coefficients de Lame
    do j = 1,NY
      do i = 1,NX-1
      value_dpressure_dx = (pressurepresent(i+1,j) - pressurepresent(i,j)) / DELTAX

      memory_dpressure_dx(i,j) = b_x_half(i) * memory_dpressure_dx(i,j) + a_x_half(i) * value_dpressure_dx

      value_dpressure_dx = value_dpressure_dx / K_x_half(i) + memory_dpressure_dx(i,j)

      pressure_xx(i,j) = value_dpressure_dx / rho(i,j)
      enddo
    enddo

    do j = 1,NY-1
      do i = 1,NX
      value_dpressure_dy = (pressurepresent(i,j+1) - pressurepresent(i,j)) / DELTAY

      memory_dpressure_dy(i,j) = b_y_half(j) * memory_dpressure_dy(i,j) + a_y_half(j) * value_dpressure_dy

      value_dpressure_dy = value_dpressure_dy / K_y_half(j) + memory_dpressure_dy(i,j)

      pressure_yy(i,j) = value_dpressure_dy / rho(i,j)
      enddo
    enddo

! calculer les deuxiemes derivees spatiales

! pour la mise a jour de pressure ci-dessous
    do j = 1,NY
      do i = 2,NX
      value_dpressurexx_dx = (pressure_xx(i,j) - pressure_xx(i-1,j)) / DELTAX

      memory_dpressurexx_dx(i,j) = b_x(i) * memory_dpressurexx_dx(i,j) + a_x(i) * value_dpressurexx_dx

      value_dpressurexx_dx = value_dpressurexx_dx / K_x(i) + memory_dpressurexx_dx(i,j)

      dpressurexx_dx(i,j) = value_dpressurexx_dx
      enddo
    enddo

! pour la mise a jour ci-dessous
    do j = 2,NY
      do i = 1,NX
      value_dpressureyy_dy = (pressure_yy(i,j) - pressure_yy(i,j-1)) / DELTAY

      memory_dpressureyy_dy(i,j) = b_y(j) * memory_dpressureyy_dy(i,j) + a_y(j) * value_dpressureyy_dy

      value_dpressureyy_dy = value_dpressureyy_dy / K_y(j) + memory_dpressureyy_dy(i,j)

      dpressureyy_dy(i,j) = value_dpressureyy_dy
      enddo
    enddo

! appliquer le schema d'evolution en temps
! on l'applique partout y compris sur certains points du bord qui n'ont pas ete calcules
! ci-dessus, ce qui est faux, mais ce n'est pas grave car on efface ces fausses valeurs
! juste apres en appliquant les conditions de Dirichlet ci-dessous
  pressurefuture(:,:) = - pressurepast(:,:) + 2.d0 * pressurepresent(:,:) + &
                                  DELTAT*DELTAT * (dpressurexx_dx(:,:) + dpressureyy_dy(:,:)) * kappa(:,:)

! imposer les conditions de bord rigide de Dirichlet (pression nulle)
! this applies Dirichlet at the bottom of the C-PML layers,
! which is the right condition to implement in order for C-PML to remain stable at long times

! bord de gauche
  pressurefuture(1,:) = ZERO

! bord de droite
  pressurefuture(NX,:) = ZERO

! bord du bas
  pressurefuture(:,1) = ZERO

! bord du haut
  pressurefuture(:,NY) = ZERO

! ajouter la source (Ricker) au point de la grille ou est situee la pression source
    a = pi*pi*f0*f0
    t = dble(it-1)*DELTAT
    force_x = factor * (1.d0-2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)
    force_y = factor * (1.d0-2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)
    pressurefuture(i_source,j_source) = pressurefuture(i_source,j_source) + force_x / rho(i_source,j_source)

! store seismograms
  do irec = 1,NREC
    sispressure(it,irec) = pressurepresent(ix_rec(irec),iy_rec(irec))
  enddo

! output information every IT_AFFICHE time steps, and at time it=5
  if(mod(it,IT_AFFICHE) == 0 .or. it == 5) then

! print max absolute value of pressure
    pressure_max_all = maxval(abs(pressurepresent))
    print *,'time step it, time t = ',it,dble(it-1)*DELTAT
    print *,'max absolute value of pressure is ',pressure_max_all

! check stability of the code, exit if unstable
    if(pressure_max_all > STABILITY_THRESHOLD) then
      stop 'code became unstable and blew up'
    endif

! display de la pression sous forme d'image en couleur
    call create_color_image(pressurepresent,NX,NY,it)
    print *,'image file written'
    print *

  endif

! move new values to old values (the present becomes the past, the future becomes the present)
  pressurepast(:,:) = pressurepresent(:,:)
  pressurepresent(:,:) = pressurefuture(:,:)

  enddo   ! end of the time loop

! save seismograms
  call write_seismograms(sispressure,NSTEP,NREC,DELTAT)

  end program seismic_CPML_2D_pressure


!----
!----  routine de sauvegarde des sismogrammes
!----

  subroutine write_seismograms(sispressure,nt,NREC,DELTAT)

! save the seismograms in ASCII format

  implicit none

  integer nt,NREC
  double precision DELTAT

  double precision, dimension(nt,NREC) :: sispressure

  integer irec,it

  character(len=100) nom_fichier

! pressure component
  do irec=1,NREC

    write(nom_fichier,"('pressure_file_',i3.3,'.dat')") irec
    open(unit=11,file=nom_fichier,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT),' ',sngl(sispressure(it,irec))
    enddo
    close(11)

  enddo

  end subroutine write_seismograms


!----
!----  routine d'affichage de la pression sous forme d'image en couleurs
!----

  subroutine create_color_image(donnees_image_color_2D,NX,NY,it)

! to display the snapshots : " display image*.gif " or " gimp image*.gif "

  implicit none

! threshold in percent of maximum amplitude below which we do not display anything
  double precision, parameter :: CUTVECT = 0.01d0

! display non lineaire pour rehausser les faibles amplitudes sur les images couleur
  double precision, parameter :: POWER_DISPLAY_COLOR = 0.30d0

  integer NX,NY,it

  double precision, dimension(NX,NY) :: donnees_image_color_2D

  integer ix,iy,R,G,B,dixmilliers,milliers,centaines,dizaines,unites,reste,current_rec

  double precision amplitude_max,valeur_normalisee

  character(len=100) nom_fichier

! create temporary image files in binary PNM P6 format (smaller) or ASCII PNM P3 format (easier to edit)
  logical, parameter :: BINARY_FILE = .false.

! ASCII code of character '0' and of carriage return character
  integer, parameter :: ascii_code_of_zero = 48, ascii_code_of_carriage_return = 10

! ouverture du fichier image
  write(nom_fichier,"('image',i6.6,'.pnm')") it

! first delete the file, just in case it was previously bigger
  open(unit=27,file=nom_fichier,status='unknown')
  close(unit=27,status='delete')

! ouvrir le fichier
  if(BINARY_FILE) then

    open(unit=27,file=nom_fichier,status='unknown',access='direct',recl=1)
    write(27,rec=1) 'P'
    write(27,rec=2) '6' ! ecrire P6 = format d'image PNM binaire
    write(27,rec=3) char(ascii_code_of_carriage_return)

! compute and write horizontal size
    reste = NX

    dixmilliers = reste / 10000
    reste = reste - 10000 * dixmilliers

    milliers = reste / 1000
    reste = reste - 1000 * milliers

    centaines = reste / 100
    reste = reste - 100 * centaines

    dizaines = reste / 10
    reste = reste - 10 * dizaines

    unites = reste

    if(dixmilliers == 0) then
      write(27,rec=4) ' '
    else
      write(27,rec=4) char(dixmilliers + ascii_code_of_zero)
    endif

    if(dixmilliers == 0 .and. milliers == 0) then
      write(27,rec=5) ' '
    else
      write(27,rec=5) char(milliers + ascii_code_of_zero)
    endif

    if(dixmilliers == 0 .and. milliers == 0 .and. centaines == 0) then
      write(27,rec=6) ' '
    else
      write(27,rec=6) char(centaines + ascii_code_of_zero)
    endif

    write(27,rec=7) char(dizaines + ascii_code_of_zero)
    write(27,rec=8) char(unites + ascii_code_of_zero)
    write(27,rec=9) ' '

! compute and write vertical size
    reste = NY

    dixmilliers = reste / 10000
    reste = reste - 10000 * dixmilliers

    milliers = reste / 1000
    reste = reste - 1000 * milliers

    centaines = reste / 100
    reste = reste - 100 * centaines

    dizaines = reste / 10
    reste = reste - 10 * dizaines

    unites = reste

    if(dixmilliers == 0) then
      write(27,rec=10) ' '
    else
      write(27,rec=10) char(dixmilliers + ascii_code_of_zero)
    endif

    if(dixmilliers == 0 .and. milliers == 0) then
      write(27,rec=11) ' '
    else
      write(27,rec=11) char(milliers + ascii_code_of_zero)
    endif

    if(dixmilliers == 0 .and. milliers == 0 .and. centaines == 0) then
      write(27,rec=12) ' '
    else
      write(27,rec=12) char(centaines + ascii_code_of_zero)
    endif

    write(27,rec=13) char(dizaines + ascii_code_of_zero)
    write(27,rec=14) char(unites + ascii_code_of_zero)
    write(27,rec=15) char(ascii_code_of_carriage_return)

! nombre de nuances
    write(27,rec=16) '2'
    write(27,rec=17) '5'
    write(27,rec=18) '5'
    write(27,rec=19) char(ascii_code_of_carriage_return)

! block of image data starts at sixteenth character
    current_rec = 20

  else

    open(unit=27,file=nom_fichier,status='unknown')
    write(27,"('P3')") ! ecrire P3 = format d'image PNM ASCII
    write(27,*) NX,NY  ! ecrire la taille
    write(27,*) '255'  ! nombre de nuances

  endif

! calculer l'amplitude maximum
  amplitude_max = maxval(abs(donnees_image_color_2D))

! dans le format PNM, l'image commence par le coin en haut a gauche
  do iy=NY,1,-1
    do ix=1,NX

! supprimer les petites amplitudes considerees comme du bruit
      if(abs(donnees_image_color_2D(ix,iy)) < amplitude_max * CUTVECT) then

! use black background where amplitude is negligible
          R = 0
          G = 0
          B = 0

      else

! definir les donnees comme etant le champ normalise entre [-1:1]
! et converti a l'entier le plus proche
! en se rappelant que l'amplitude peut etre negative
        valeur_normalisee = donnees_image_color_2D(ix,iy) / amplitude_max

! supprimer valeurs en dehors de [-1:+1]
        if(valeur_normalisee < -1.d0) then
          valeur_normalisee = -1.d0
        endif
        if(valeur_normalisee > 1.d0) then
          valeur_normalisee = 1.d0
        endif

! utiliser rouge si champ positif, bleu si negatif, pas de vert
        if(valeur_normalisee >= 0.d0) then
          R = nint(255.d0*valeur_normalisee**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(valeur_normalisee)**POWER_DISPLAY_COLOR)
        endif

     endif

! ecrire l'image en couleur
      if(BINARY_FILE) then

! first write red
        write(27,rec=current_rec) char(R)
        current_rec = current_rec + 1

! then write green
        write(27,rec=current_rec) char(G)
        current_rec = current_rec + 1

! then write blue
        write(27,rec=current_rec) char(B)
        current_rec = current_rec + 1

      else

        write(27,*) R,G,B

      endif

    enddo
  enddo

! fermer le fichier
  close(27)

  end subroutine create_color_image

