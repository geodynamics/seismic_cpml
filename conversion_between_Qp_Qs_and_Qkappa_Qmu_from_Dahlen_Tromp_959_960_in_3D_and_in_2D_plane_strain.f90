
  program conversion

! Dimitri Komatitsch, CNRS Marseille, France, July 2018

! see formulas 9.59 and 9.60 in the book of Dahlen and Tromp, 1998
! (in that book, P is called alpha and S is called beta).
! See also file formulas_to_convert_between_Qkappa_Qmu_and_Qp_Qs_in_3D_and_in_2D_plane_strain.pdf in this directory.

  implicit none

  integer :: iconversion_type
  integer :: idimension

  double precision :: Qkappa,Qmu,Qp,Qs,cp,cs
  double precision :: inverse_of_Qp,inverse_of_Qmu,inverse_of_Qkappa,COEFFICIENT

  print *,'1 = you want to perform the conversion in 3D'
  print *,'2 = you want to perform the conversion in 2D plane strain'
  read(*,*) idimension
  if (idimension < 1 .or. idimension > 2) stop 'error: incorrect value of idimension'
  if (idimension == 1) then
    COEFFICIENT = 4.d0 / 3.d0
  else
    COEFFICIENT = 1.d0
  endif
  print *

  print *,'1 = you want to convert from (Qp,Qs) to (QKappa,Qmu)'
  print *,'2 = you want to convert from (QKappa,Qmu) to (Qp,Qs)'
  read(*,*) iconversion_type
  if (iconversion_type < 1 .or. iconversion_type > 2) stop 'error: incorrect value of iconversion_type'
  print *

! get the input values from the user
  if (iconversion_type == 1) then
    print *,'enter Qp:'
    read(*,*) Qp
    print *,'enter Qs:'
    read(*,*) Qs
    print *
  else
    print *,'enter QKappa:'
    read(*,*) QKappa
    print *,'enter Qmu:'
    read(*,*) Qmu
    print *
  endif

! enter the cp and cs velocities of the medium, at the frequency at which you want this conversion to be performed
  print *,'enter the cp and cs velocities of the medium, at the frequency at which you want this conversion to be performed:'
  print *,'enter cp:'
  read(*,*) cp
  print *,'enter cs:'
  read(*,*) cs
  print *

  if (iconversion_type == 1) then

! Qmu is always the same as Qs
    Qmu = Qs

! for QKappa the formula is more complex
    inverse_of_Qp = 1.d0 / Qp
    inverse_of_Qmu = 1.d0 / Qmu

    inverse_of_Qkappa = (inverse_of_Qp - COEFFICIENT*(cs**2)/(cp**2) * inverse_of_Qmu) / (1.d0 - COEFFICIENT*(cs**2)/(cp**2))

    Qkappa = 1.d0/inverse_of_Qkappa

! print the result
    print *,'Qkappa = ',Qkappa
    print *,'Qmu = ',Qmu

  else ! if (iconversion_type == 2) then

! Qs is always the same as Qmu
    Qs = Qmu

! for Qp the formula is more complex
    inverse_of_Qp = (1.d0 - COEFFICIENT*(cs**2)/(cp**2))/Qkappa + COEFFICIENT*(cs**2)/(cp**2)/Qmu
    Qp = 1.d0/inverse_of_Qp

! print the result
    print *,'Qp = ',Qp
    print *,'Qs = ',Qs

  endif

  end program conversion

