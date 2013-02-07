#
# Makefile for SEISMIC_CPML Version 1.1.1, November 2009.
# Dimitri Komatitsch
# Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France
# 
SHELL=/bin/sh

O = obj

# the MEDIUM_MEMORY flag is for large 3D runs, which need more than 2 GB of memory

# Portland
#F90 = pgf90
#MPIF90 = mpif90
#FLAGS = -fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -fastsse -tp amd64e -Msmart
#MEDIUM_MEMORY = -mcmodel=medium
#OPEN_MP = -mp

# Intel (leave option -ftz, which can be *critical* for performance)
#F90 = ifort
#MPIF90 = mpif90
#FLAGS = -O3 -xHost -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -ftz
#MEDIUM_MEMORY = -mcmodel=medium
#OPEN_MP = -openmp -openmp-report1

# IBM xlf
#F90 = xlf_r
#MPIF90 = mpxlf_r
#FLAGS = -O3 -qfree=f90 -qhalt=w -qsave
#MEDIUM_MEMORY = -q64
#OPEN_MP = -qsmp=omp

# GNU gfortran
F90 = gfortran
MPIF90 = mpif90
FLAGS = -std=f2003 -fimplicit-none -frange-check -O3 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow
MEDIUM_MEMORY = -mcmodel=medium
#OPEN_MP = -fopenmp

default: clean seismic_CPML_2D_isotropic_second_order seismic_CPML_2D_isotropic_fourth_order seismic_CPML_2D_anisotropic seismic_PML_Collino_2D_isotropic seismic_PML_Collino_3D_isotropic_OpenMP seismic_CPML_3D_isotropic_MPI_OpenMP seismic_CPML_2D_poroelastic_fourth_order seismic_CPML_3D_viscoelastic_MPI seismic_PML_Collino_2D_anisotropic_fourth seismic_ADEPML_2D_elastic_RK4_eighth_order seismic_ADEPML_2D_viscoelastic_RK4_eighth_order

all: default

clean:
	/bin/rm -f *.o xseismic_CPML_2D_isotropic_second_order xseismic_CPML_2D_isotropic_fourth_order xseismic_CPML_2D_anisotropic xseismic_PML_Collino_2D_isotropic xseismic_CPML_3D_isotropic_MPI_OpenMP xseismic_PML_Collino_3D_isotropic_OpenMP xseismic_CPML_2D_poroelastic_fourth_order xseismic_CPML_3D_viscoelastic_MPI xseismic_PML_Collino_2D_anisotropic_fourth xseismic_ADEPML_2D_elastic_RK4_eighth_order xseismic_ADEPML_2D_viscoelastic_RK4_eighth_order

seismic_ADEPML_2D_elastic_RK4_eighth_order:
	$(F90) $(FLAGS) -o xseismic_ADEPML_2D_elastic_RK4_eighth_order seismic_ADEPML_2D_elastic_RK4_eighth_order.f90

seismic_ADEPML_2D_viscoelastic_RK4_eighth_order:
	$(F90) $(FLAGS) -o xseismic_ADEPML_2D_viscoelastic_RK4_eighth_order seismic_ADEPML_2D_viscoelastic_RK4_eighth_order.f90

seismic_CPML_2D_poroelastic_fourth_order:
	$(F90) $(FLAGS) -o xseismic_CPML_2D_poroelastic_fourth_order seismic_CPML_2D_poroelastic_fourth_order.f90

seismic_CPML_2D_isotropic_second_order:
	$(F90) $(FLAGS) -o xseismic_CPML_2D_isotropic_second_order seismic_CPML_2D_isotropic_second_order.f90

seismic_CPML_2D_isotropic_fourth_order:
	$(F90) $(FLAGS) -o xseismic_CPML_2D_isotropic_fourth_order seismic_CPML_2D_isotropic_fourth_order.f90

seismic_CPML_2D_anisotropic:
	$(F90) $(FLAGS) -o xseismic_CPML_2D_anisotropic seismic_CPML_2D_anisotropic.f90

seismic_PML_Collino_2D_isotropic:
	$(F90) $(FLAGS) -o xseismic_PML_Collino_2D_isotropic seismic_PML_Collino_2D_isotropic.f90

seismic_PML_Collino_2D_anisotropic_fourth:
	$(F90) $(FLAGS) -o xseismic_PML_Collino_2D_anisotropic_fourth seismic_PML_Collino_2D_anisotropic_fourth.f90

seismic_PML_Collino_3D_isotropic_OpenMP:
	$(F90) $(FLAGS) $(MEDIUM_MEMORY) $(OPEN_MP) -o xseismic_PML_Collino_3D_isotropic_OpenMP seismic_PML_Collino_3D_isotropic_OpenMP.f90

seismic_CPML_3D_isotropic_MPI_OpenMP:
	$(MPIF90) $(FLAGS) $(MEDIUM_MEMORY) $(OPEN_MP) -o xseismic_CPML_3D_isotropic_MPI_OpenMP seismic_CPML_3D_isotropic_MPI_OpenMP.f90

seismic_CPML_3D_viscoelastic_MPI:
	$(MPIF90) $(FLAGS) $(MEDIUM_MEMORY) $(OPEN_MP) -o xseismic_CPML_3D_viscoelastic_MPI seismic_CPML_3D_viscoelastic_MPI.f90

