PDE_Parameters.diameter = 1032*1e-9;                            %Scaling parameter for the nanoparticles
PDE_Parameters.particle_Density = 1009;                         %Density of particles (kg/m^3)
PDE_Parameters.affinity = 2.97e-9/PDE_Parameters.dx;            %Boundary condition flux condition parameter
PDE_Parameters.K = 60;                                          %Carrying capacity NP/Cell
PDE_Parameters.num_Cells = 100000;
PDE_Parameters.cell_Area = 830e-12;
PDE_Parameters.initial_Particle_Count = 1e7;                                %Number of particles at the beginning of the experiment