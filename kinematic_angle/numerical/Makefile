F90    = gfortran
FFLAGS = -Wall -g -ffpe-trap=invalid,zero
PDFLIB = `cernlib pdflib`

all: fig1.exe grv_figures.exe main_deflection_by.exe main_deflection_block.exe\
 main_deflection_bodek.exe

bin/%.o: %.f90
	${F90} -J bin ${FFLAGS} $< -c -o $@

fig1.exe: $(addprefix bin/,fig1.o by03.o grv94.o ogpf.o)
	${F90} ${FFLAGS} $^ ${PDFLIB} -o $@
grv_figures.exe: $(addprefix bin/,grv_figures.o grv94.o ogpf.o)
	${F90} ${FFLAGS} $^ ${PDFLIB} -o $@
main_deflection_by.exe: $(addprefix bin/,main_deflection_by.o grv94.o by03.o\
 multidim_integrate.o kinematic.o)
	${F90} ${FFLAGS} $^ ${PDFLIB} -o $@
main_deflection_bodek.exe: $(addprefix bin/,main_deflection_bodek.o grv98.o\
 bodek2021.o multidim_integrate.o kinematic.o ogpf.o)
	${F90} ${FFLAGS} $^ ${PDFLIB} -o $@
main_deflection_block.exe: $(addprefix bin/,main_deflection_block.o block.o\
 multidim_integrate.o kinematic.o)
	${F90} ${FFLAGS} $^ ${PDFLIB} -o $@

-include dep.mk

dep.mk: fig1.f90 by03.f90 grv94.f90 ogpf.f90 main_deflection_block.f90\
 main_deflection_by.f90 multidim_integrate.f90 kinematic.f90 block.f90\
 main_deflection_bodek.f90 grv98.f90 bodek2021.f90
	makedepf90 -b bin $^ > $@

clean:
	rm -rf dep.mk bin/* *.exe

.PHONY: all clean
