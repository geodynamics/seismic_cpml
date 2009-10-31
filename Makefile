#
# Makefile
#
# Dimitri Komatitsch, University of Pau, April 2007
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

# Intel
#F90 = ifort
#MPIF90 = mpif90
#FLAGS = -O3 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -assume byterecl -check nobounds
#MEDIUM_MEMORY = -mcmodel=medium -i-dynamic
#OPEN_MP = -openmp

# IBM xlf
#F90 = xlf_r
#MPIF90 = mpxlf_r
#FLAGS = -O3 -qfree=f90 -qhalt=w -qsave
#MEDIUM_MEMORY =
#OPEN_MP = -qsmp=omp

# GNU gfortran
F90 = gfortran
MPIF90 = mpif90
FLAGS = -std=gnu -fimplicit-none -frange-check -O2 -Wunused-labels -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow
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

