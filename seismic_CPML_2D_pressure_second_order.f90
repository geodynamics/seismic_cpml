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
