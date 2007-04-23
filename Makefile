#
# Makefile
#
# Dimitri Komatitsch, University of Pau, April 2007
# 
SHELL=/bin/sh

O = obj

# Portland
#F90 = pgf90
#MPIF90 = mpif90
#FLAGS=-fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -fastsse -tp amd64e -Msmart
#OPEN_MP=-mp

# Intel
#F90 = ifort
#MPIF90 = mpif90
#FLAGS=-O3 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -assume byterecl -check nobounds
#OPEN_MP=-openmp

# GNU gfortran
F90 = gfortran
MPIF90 = mpif90
FLAGS = -std=gnu -fimplicit-none -frange-check -O2 -Wunused-labels -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow
OPEN_MP=-fopenmp

# large 3D runs need more than 2 GB of memory
MEDIUM_MEMORY=-mcmodel=medium

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

