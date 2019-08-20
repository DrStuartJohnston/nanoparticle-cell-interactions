%% Script to define parameters used in the traditional dosage model.

addpath('./Particle Properties');

particle_Type = '282nm_Coreshell_Hela';                                 %Choice of nanoparticle-cell pair

PDE_Parameters.num_Nodes = 10001;                                       %Number of spatial nodes
PDE_Parameters.width = 2.62*1e-3;                                       %Depth of media
PDE_Parameters.dx = PDE_Parameters.width/(PDE_Parameters.num_Nodes-1);  %Grid spacing for numerical solution
PDE_Parameters.dt = 1;                                                  %Temporal spacing for numerical solution

%% Define nanoparticle-cell pair properties
if strcmpi(particle_Type,'150nm_Coreshell_Raw')
    Parameters_150nm_Coreshell_Raw;
elseif strcmpi(particle_Type,'282nm_Coreshell_Hela')
    Parameters_282nm_Coreshell_Hela;
elseif strcmpi(particle_Type,'1032nm_Capsule_Raw')
    Parameters_1032nm_Capsule_Raw;
else
    error('Incorrect Particle Name')
end

PDE_Parameters.avogrado_Number = 6.022*1e23;                        %Avogrado's number
PDE_Parameters.gas_Constant = 8.3144598;                            %Gas constant
PDE_Parameters.media_Density = 1000;                                %Density of media that particles are immersed in
PDE_Parameters.media_Viscosity = 0.00101;                           %Viscosity of media that particles are immersed in
PDE_Parameters.temperature = 310.15;                                %Temperature of fluid in K
PDE_Parameters.gravity = 9.81;                                      %Gravitational acceleration
PDE_Parameters.final_Time = 24*3600;                                %Final time of experiment
PDE_Parameters.r = PDE_Parameters.affinity*PDE_Parameters.dx;       %Affinity Parameter
PDE_Parameters.dish_Area = pi*(0.0221/2)^2;                         %Area of base of culture well
PDE_Parameters.S = PDE_Parameters.num_Cells*PDE_Parameters.cell_Area/PDE_Parameters.dish_Area;  %Level of confluency                                             %Surface area parameter
PDE_Parameters.volume = PDE_Parameters.dish_Area*PDE_Parameters.width;                          %Volume of culture well
PDE_Parameters.particle_Concentration = PDE_Parameters.initial_Particle_Count/PDE_Parameters.volume;    %Initial nanoparticle concentration
PDE_Parameters.carrying_Capacity_Ratio = PDE_Parameters.initial_Particle_Count/PDE_Parameters.K;        %Number of nanoparticles relative to the total carrying capacity          
nondim = 'yes';

PDE_Parameters.D = PDE_Parameters.gas_Constant*PDE_Parameters.temperature./(3*PDE_Parameters.avogrado_Number*pi*PDE_Parameters.media_Viscosity*PDE_Parameters.diameter);    %Diffusion coefficient
PDE_Parameters.V = PDE_Parameters.gravity*(PDE_Parameters.particle_Density-PDE_Parameters.media_Density).*PDE_Parameters.diameter.^2/(18*PDE_Parameters.media_Viscosity);   %Sedimentation coefficient