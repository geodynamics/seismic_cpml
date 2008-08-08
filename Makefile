#
# Makefile
#
# Version 1.0
# Dimitri Komatitsch
# Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France, April 2007
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

# Intel (leave option -ftz, which is *critical* for performance)
#F90 = ifort
#MPIF90 = mpif90
#FLAGS = -O3 -xP -vec-report0 -e03 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -fpe3 -ftz
#MEDIUM_MEMORY = -mcmodel=medium
#OPEN_MP = -openmp -openmp-report1

# IBM xlf
#F90 = xlf_r
#MPIF90 = mpxlf_r
#FLAGS = -O3 -qfree=f90 -qhalt=w -qsave
#MEDIUM_MEMORY =
#OPEN_MP = -qsmp=omp

# GNU gfortran
F90 = gfortran
MPIF90 = mpif90
FLAGS = -std=f2003 -fimplicit-none -frange-check -O3 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fno-trapping-math
MEDIUM_MEMORY = -mcmodel=medium
OPEN_MP = -fopenmp

default: clean seismic_CPML_2D_iso seismic_CPML_2D_aniso seismic_PML_Collino_2D_iso seismic_PML_Collino_3D_iso_OpenMP seismic_CPML_3D_iso_MPI_OpenMP

all: default

clean:
	/bin/rm -f *.o xseismic_CPML_2D_iso xseismic_CPML_2D_aniso xseismic_PML_Collino_2D_iso xseismic_CPML_3D_iso_MPI_OpenMP xseismic_PML_Collino_3D_iso_OpenMP

seismic_CPML_2D_iso:
	$(F90) $(FLAGS) -o xseismic_CPML_2D_iso seismic_CPML_2D_iso.f90

seismic_CPML_2D_aniso:
	$(F90) $(FLAGS) -o xseismic_CPML_2D_aniso seismic_CPML_2D_aniso.f90

seismic_PML_Collino_2D_iso:
	$(F90) $(FLAGS) -o xseismic_PML_Collino_2D_iso seismic_PML_Collino_2D_iso.f90

seismic_PML_Collino_3D_iso_OpenMP:
	$(F90) $(FLAGS) $(MEDIUM_MEMORY) $(OPEN_MP) -o xseismic_PML_Collino_3D_iso_OpenMP seismic_PML_Collino_3D_iso_OpenMP.f90

seismic_CPML_3D_iso_MPI_OpenMP:
	$(MPIF90) $(FLAGS) $(MEDIUM_MEMORY) $(OPEN_MP) -o xseismic_CPML_3D_iso_MPI_OpenMP seismic_CPML_3D_iso_MPI_OpenMP.f90

