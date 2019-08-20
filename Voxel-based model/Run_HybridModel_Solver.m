%% Define relevant inputs and calculate numerical solution for hybrid model

rng(1); %Fix rng seed for consistent cell characteristic distribution sampling

x = 0:PDE_Parameters.dx:PDE_Parameters.width;                                                                   %Spatial domain

initial_Concentration_Value = 1.0;                                                                  %Value of initial condition
initial_Concentration_Width = 1.0*PDE_Parameters.width;                                             %Width of initial condition

initial_Concentration = zeros(1,PDE_Parameters.num_Nodes);                                          %Initialise initial condition vector

for i = 1:PDE_Parameters.num_Nodes
    if x(i) <= initial_Concentration_Width
        initial_Concentration(i) = initial_Concentration_Value;                                     %Define initial condition
    end
end

%% Nondimensionalisation

if strcmpi(nondim,'yes')
    PDE_Parameters.rescaled_x = x/PDE_Parameters.width;                                                         %Rescaled domain
    PDE_Parameters.rescaled_Final_Time = PDE_Parameters.final_Time*abs(PDE_Parameters.V)/PDE_Parameters.width;  %Rescaled end of experiment time
    PDE_Parameters.rescaled_dt = PDE_Parameters.dt*abs(PDE_Parameters.V)/PDE_Parameters.width;                  %Rescaled time step
    PDE_Parameters.rescaled_D = PDE_Parameters.D/(abs(PDE_Parameters.V)*PDE_Parameters.width);                  %Rescaled diffusion coefficient
    PDE_Parameters.rescaled_V = 1;                                                                              %Rescaled sedimentation coefficient
    PDE_Parameters.rescaled_dx = PDE_Parameters.dx/PDE_Parameters.width;                                        %Rescaled grid spacing
    PDE_Parameters.rescaled_S = PDE_Parameters.S/PDE_Parameters.V;
else
    PDE_Parameters.rescaled_x = x;                                                       %Rescaled domain
    PDE_Parameters.rescaled_Final_Time = PDE_Parameters.final_Time;                      %Rescaled end of experiment time
    PDE_Parameters.rescaled_dt = PDE_Parameters.dt;                                      %Rescaled time step
    PDE_Parameters.rescaled_D = PDE_Parameters.D;                                        %Rescaled diffusion coefficient
    PDE_Parameters.rescaled_V = PDE_Parameters.V;                                        %Rescaled sedimentation coefficient
    PDE_Parameters.rescaled_dx = PDE_Parameters.dx;                                      %Rescaled grid spacing
    PDE_Parameters.rescaled_S = PDE_Parameters.S;
    initial_Concentration = initial_Concentration*PDE_Parameters.particle_Concentration; %Scale initial concentration
end

PDE_Parameters.initial = [initial_Concentration,0];                                      %Pass in initial condition

time = 0:PDE_Parameters.rescaled_dt:PDE_Parameters.rescaled_Final_Time;                  %Vector of solution time points

C = initial_Concentration';                                                              %Initial solution
total_Associated = zeros(numel(time),1);                                                 %Vector of number of associated nanoparticles 
[aVector,bVector,cVector] = define_ThomasAlg_Vectors(PDE_Parameters);                    %Calculate vectors for Thomas' Algorithm 
opts = [];                                            

cell_Locations = randperm(IBM_Parameters.total_Cell_Spaces,IBM_Parameters.number_of_Cells); %Define cell locations

%% Generate nanoparticle-cell affinity distribution if needed
if strcmpi(IBM_Parameters.affinity_Variable,'yes')
    cell_Affinity = lognrnd(log(PDE_Parameters.r)-IBM_Parameters.affinity_StDev^2/2,IBM_Parameters.affinity_StDev,IBM_Parameters.total_Cell_Spaces,1);
else  
    cell_Affinity = PDE_Parameters.r;
end

%% Generate cell carrying capacity distribution if needed
if strcmpi(IBM_Parameters.capacity_Variable,'yes')
    cell_Capacity = IBM_Parameters.cell_Capacity*lognrnd(log(1)-IBM_Parameters.capacity_StDev^2/2,IBM_Parameters.capacity_StDev,IBM_Parameters.total_Cell_Spaces,1);
    max_NP = sum(cell_Capacity(cell_Locations));    %Maximum number of associated nanoparticles
else
    cell_Capacity = IBM_Parameters.cell_Capacity;
    max_NP = IBM_Parameters.number_of_Cells*cell_Capacity; %Maximum number of associated nanoparticles
end

%% Generate correlated affinity and capacity distributions if needed
if strcmpi(IBM_Parameters.related_Capacity_Affinity,'yes')
    tmp_Distribution = lognrnd(log(1)-IBM_Parameters.capacity_StDev^2/2,IBM_Parameters.capacity_StDev,IBM_Parameters.total_Cell_Spaces,1);  %Temporary distribution
    cell_Capacity = IBM_Parameters.cell_Capacity*tmp_Distribution;                                                          %Capacity distribution
    tmp_CDF = logncdf(tmp_Distribution,log(1)-IBM_Parameters.capacity_StDev^2/2,IBM_Parameters.capacity_StDev);             %Corresponding CDF
    shifted_Distribution = logninv(tmp_CDF,log(1)-IBM_Parameters.affinity_StDev^2/2,IBM_Parameters.affinity_StDev);         %Distribution values for affinity distribution
    cell_Affinity = PDE_Parameters.r*shifted_Distribution;                                                                  %Affinity distribution
    max_NP = sum(cell_Capacity(cell_Locations));                                                                            %Maximum number of associated nanoparticles
end

NPs_per_Cell = zeros(IBM_Parameters.number_of_Cells,1);         %Vector of nanoparticles per cell

if strcmpi(IBM_Parameters.capacity_Variable,'yes') && strcmpi(IBM_Parameters.affinity_Variable,'yes')
    dCell = cell_Affinity(cell_Locations).*((cell_Capacity(cell_Locations)-NPs_per_Cell)./cell_Capacity(cell_Locations));   %Flux through the boundary for each cell
elseif strcmpi(IBM_Parameters.capacity_Variable,'yes') && strcmpi(IBM_Parameters.affinity_Variable,'no')
    dCell = cell_Affinity.*((cell_Capacity(cell_Locations)-NPs_per_Cell)./cell_Capacity(cell_Locations)); %Flux through the boundary for each cell
elseif strcmpi(IBM_Parameters.capacity_Variable,'no') && strcmpi(IBM_Parameters.affinity_Variable,'yes')
    dCell = cell_Affinity(cell_Locations).*((cell_Capacity-NPs_per_Cell)./cell_Capacity); %Flux through the boundary for each cell
end
mean_dCell = mean(dCell); %Mean flux through the boundary 

heterogeneity_Evolution = zeros(floor(numel(time)/1000),1); %Evolution of the apparent heterogeneity
NP_per_Cell_Evolution = zeros(floor(numel(time)/1000),IBM_Parameters.number_of_Cells);  %Evolution of the number of nanoparticles per cell

for tCounter = 1:numel(time)
    if mod(tCounter,1000) == 0
        tCounter/numel(time)
        tmp_Parameters = fitdist(NPs_per_Cell,'lognormal');                     %Fit distribution to number of nanoparticles per cell
        heterogeneity_Evolution(floor(tCounter/1000)) = tmp_Parameters.sigma;   %Evolution of apparent heterogeneity
        NP_per_Cell_Evolution(floor(tCounter/1000),:) = NPs_per_Cell;           %Evolution of the number of nanoparticles per cell
    end
    total_Associated(tCounter) = (sum(PDE_Parameters.initial)-sum(C))/sum(PDE_Parameters.initial)/PDE_Parameters.num_Cells*PDE_Parameters.initial_Particle_Count; %Number of associated nanoparticles
    bVector(end) = update_bVector(PDE_Parameters,total_Associated(tCounter),mean_dCell);    %Update vector for Thomas' Algorithm corresponding to boundary condition
    C = Thomas_Algorithm(aVector,bVector,cVector,C);                                        %Solve tridiagonal system via Thomas' algorithm
    if tCounter > 1
        NPs_per_Cell = NPs_per_Cell + (total_Associated(tCounter)-total_Associated(tCounter-1))*IBM_Parameters.number_of_Cells*dCell/sum(dCell); %Split the flux through the boundary into individual cells
    else
        NPs_per_Cell = NPs_per_Cell + total_Associated(tCounter)*IBM_Parameters.number_of_Cells*dCell/sum(dCell);                                %Split the flux through the boundary into individual cells
    end
    if strcmpi(IBM_Parameters.capacity_Variable,'yes') && strcmpi(IBM_Parameters.affinity_Variable,'yes')
        dCell = cell_Affinity(cell_Locations).*((cell_Capacity(cell_Locations)-NPs_per_Cell)./cell_Capacity(cell_Locations)); %Flux through the boundary for each cell
    elseif strcmpi(IBM_Parameters.capacity_Variable,'yes') && strcmpi(IBM_Parameters.affinity_Variable,'no')
        dCell = cell_Affinity.*((cell_Capacity(cell_Locations)-NPs_per_Cell)./cell_Capacity(cell_Locations)); %Flux through the boundary for each cell
    elseif strcmpi(IBM_Parameters.capacity_Variable,'no') && strcmpi(IBM_Parameters.affinity_Variable,'yes')
        dCell = cell_Affinity(cell_Locations).*((cell_Capacity-NPs_per_Cell)./cell_Capacity); %Flux through the boundary for each cell
    end
    mean_dCell = mean(dCell); %Mean flux through the boundary 

end

final_Profile = C;                      %Final solution
scaled_Final_Profile = final_Profile;   %Final solution
uptake_Curve = total_Associated;        %Nanoparticle uptake curve 
unscaled_Time = linspace(0,PDE_Parameters.final_Time/3600,numel(time)); %Unscaled time vector

figure; plot(unscaled_Time,uptake_Curve,'linewidth',2) %Plot nanoparticle uptake

figure; plot(unscaled_Time(1000:1000:end),heterogeneity_Evolution,'linewidth',2) %Plot evolution of apparent heterogeneity