%% Script to define parameters used in the voxel-based model.

IBM_Parameters.scaling = 100;                                                                                                       %Scaling parameter for efficient time-step sampling
IBM_Parameters.NPs_per_Cell_Initial = PDE_Parameters.initial_Particle_Count/PDE_Parameters.num_Cells*IBM_Parameters.scaling;        %Total number of nanoparticles relative to number of cells
IBM_Parameters.cell_Capacity = PDE_Parameters.K;                %Carrying capacity of the cells
IBM_Parameters.cell_Spaces_x = 20;                              %Number of voxels (x)
IBM_Parameters.cell_Spaces_y = 20;                              %Number of voxels (y)
IBM_Parameters.total_Cell_Spaces = IBM_Parameters.cell_Spaces_x*IBM_Parameters.cell_Spaces_y; %Number of boundary voxels
IBM_Parameters.cell_Density = PDE_Parameters.S;                                               %Level of confluency in the experiment
IBM_Parameters.number_of_Cells = ceil(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.cell_Density);    %Number of voxels containing cells
IBM_Parameters.grid_Size = sqrt(PDE_Parameters.cell_Area);                                              %Voxel size
IBM_Parameters.grid_Height = round(PDE_Parameters.width/IBM_Parameters.grid_Size);                      %Number of voxels (z)
IBM_Parameters.t_End = PDE_Parameters.final_Time;                                                       %Final time of experiment
IBM_Parameters.D = PDE_Parameters.D;                                                                    %Diffusion coefficient
IBM_Parameters.V = PDE_Parameters.V;                                                                    %Sedimentation coefficient
IBM_Parameters.affinity_Scale = PDE_Parameters.r/(PDE_Parameters.D);                                    %IBM affinity parameter
IBM_Parameters.affinity = IBM_Parameters.affinity_Scale*IBM_Parameters.grid_Size;                       %Scaled affinity parameter
IBM_Parameters.time_Step = 10;                                                                          %Tau-leaping time step
IBM_Parameters.NP_Initial_Density = ceil(IBM_Parameters.NPs_per_Cell_Initial*IBM_Parameters.number_of_Cells/(IBM_Parameters.total_Cell_Spaces*IBM_Parameters.grid_Height)); %Number of nanoparticles per voxel
IBM_Parameters.affinity_Variable = 'yes';                       %Is the affinity a constant or does it follow a distribution ('yes' or 'no')
IBM_Parameters.affinity_StDev = 0.65551;                        %Standard devation of the affinity distribution
IBM_Parameters.capacity_Variable = 'yes';                       %Is the carrying capacity a constant or does it follow a distribution ('yes' or 'no')
IBM_Parameters.capacity_StDev = 0.65551;                        %Standard deviation of the capacity distribution
IBM_Parameters.related_Capacity_Affinity = 'yes';               %Are the affinity and carrying capacity correlated ('yes' or 'no')

IBM_Parameters.distribution_Fits = ["poisson"];                 %Distribution to fit to model output ("normal";"lognormal";"binomial";"negbin")
IBM_Parameters.cell_Cycle = 'no';                               %Include cell cycle ('yes' or 'no')
IBM_Parameters.cell_Cycle_Dependent_Affinity = 'no';            %Does the affinity depend on cell cycle phase ('yes' or 'no')
if strcmpi(particle_Type,'150nm_Coreshell_Raw') || strcmpi(particle_Type,'1032nm_Capsule_Raw')
    IBM_Parameters.G1_Time = 0.55*15*60*60;             %G1 phase length
    IBM_Parameters.S_Time = 0.35*15*60*60;              %S phase length
    IBM_Parameters.G2_Time = 0.075*15*60*60;            %G2 phase length
    IBM_Parameters.M_Time = 0.025*15*60*60;             %M phase length
elseif strcmpi(particle_Type,'282nm_Coreshell_Hela')
    IBM_Parameters.G1_Time = 400*60;                    %G1 phase length
    IBM_Parameters.S_Time = 500*60;                     %S phase length
    IBM_Parameters.G2_Time = 175*60;                    %G2 phase length
    IBM_Parameters.M_Time = 50*60;                      %M phase length
end

if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
    IBM_Parameters.G1_Affinity = 0.7*IBM_Parameters.affinity;   %G1 phase affinity
    IBM_Parameters.S_Affinity = 1*IBM_Parameters.affinity;      %S phase affinity
    IBM_Parameters.G2_Affinity = 1.6*IBM_Parameters.affinity;   %G2 phase affinity
    IBM_Parameters.M_Affinity = 2*IBM_Parameters.affinity;      %M phase affinity
end

IBM_Parameters.opts = optimset('Display','off');                %Optimisation routine options - removal of display