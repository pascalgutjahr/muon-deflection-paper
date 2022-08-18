bin/fig1.o : fig1.f90 bin/ogpf.o bin/by03.o 
bin/by03.o : by03.f90 bin/grv94.o 
bin/grv94.o : grv94.f90 
bin/ogpf.o : ogpf.f90 
bin/main_deflection_block.o : main_deflection_block.f90 bin/block.o bin/kinematic.o bin/multidim_integrate.o 
bin/main_deflection_by.o : main_deflection_by.f90 bin/by03.o bin/kinematic.o bin/multidim_integrate.o 
bin/multidim_integrate.o : multidim_integrate.f90 
bin/kinematic.o : kinematic.f90 
bin/block.o : block.f90 
bin/main_deflection_bodek.o : main_deflection_bodek.f90 bin/ogpf.o bin/bodek2021.o bin/kinematic.o bin/multidim_integrate.o 
bin/grv98.o : grv98.f90 
bin/bodek2021.o : bodek2021.f90 bin/grv98.o 
