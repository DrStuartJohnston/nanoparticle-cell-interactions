PDE_Parameters.diameter = 250*1e-9;                             %Scaling parameter for the nanoparticles
PDE_Parameters.particle_Density = 1650;                        %Density of particles (kg/m^3)
PDE_Parameters.affinity = 2.14e-10/PDE_Parameters.dx;            %Boundary condition flux condition parameter
PDE_Parameters.K = 150;                                          %Carrying capacity NP/Cell
PDE_Parameters.num_Cells = 80000;
PDE_Parameters.cell_Area = 1600e-12;
PDE_Parameters.initial_Particle_Count = 4e8;                                %Number of particles at the beginning of the experiment