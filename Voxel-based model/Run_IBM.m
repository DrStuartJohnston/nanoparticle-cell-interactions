%% Run the voxel-based model.

cell_Locations = randperm(IBM_Parameters.total_Cell_Spaces,IBM_Parameters.number_of_Cells); %Define cell locations
cell_Layer = zeros(IBM_Parameters.total_Cell_Spaces,1);                                     %Define boundary layer of voxels    
cell_Layer(cell_Locations) = 1;                                                             %Define cell locations in boundary layer
NPs_per_Cell = zeros(IBM_Parameters.total_Cell_Spaces,1);                                   %Initialise number of nanoparticles per cell

NP_Domain = IBM_Parameters.NP_Initial_Density*ones(IBM_Parameters.grid_Height*IBM_Parameters.total_Cell_Spaces,1);  %Define domain of voxel-based model

extra_NPs = sum(NP_Domain) - IBM_Parameters.number_of_Cells*IBM_Parameters.NPs_per_Cell_Initial; %Ensure number of nanoparticles in the system is consistent

%% Remove extra NPs if necessary
if extra_NPs > 0    
    removed_NPs = randperm(numel(NP_Domain),extra_NPs);
    NP_Domain(removed_NPs) = NP_Domain(removed_NPs)-1;
    IBM_Parameters.NP_Initial_Density = sum(NP_Domain)/numel(NP_Domain);
end

t = 0;                              %Current time
counter = 0;                        %Define counting variable
dt = IBM_Parameters.time_Step;      %Time step in voxel-based model

%% Define nanoparticle-cell affinity distribution
if strcmpi(IBM_Parameters.affinity_Variable,'yes')
    cell_Affinity = lognrnd(log(IBM_Parameters.affinity)-IBM_Parameters.affinity_StDev^2/2,IBM_Parameters.affinity_StDev,IBM_Parameters.total_Cell_Spaces,1);
else
    cell_Affinity = IBM_Parameters.affinity;
end

%% Define cell carrying capacity distribution
if strcmpi(IBM_Parameters.capacity_Variable,'yes')
    cell_Capacity = IBM_Parameters.cell_Capacity*lognrnd(log(1)-IBM_Parameters.capacity_StDev^2/2,IBM_Parameters.capacity_StDev,IBM_Parameters.total_Cell_Spaces,1);
else
    cell_Capacity = IBM_Parameters.cell_Capacity;
end

%% Define correlated affinity and cell carrying capacity distributions
if strcmpi(IBM_Parameters.related_Capacity_Affinity,'yes')
    tmp_Distribution = lognrnd(log(1)-IBM_Parameters.capacity_StDev^2/2,IBM_Parameters.capacity_StDev,IBM_Parameters.total_Cell_Spaces,1);
    cell_Capacity = IBM_Parameters.cell_Capacity*tmp_Distribution;
    tmp_CDF = logncdf(tmp_Distribution,log(1)-IBM_Parameters.capacity_StDev^2/2,IBM_Parameters.capacity_StDev);
    shifted_Distribution = logninv(tmp_CDF,log(1)-IBM_Parameters.affinity_StDev^2/2,IBM_Parameters.affinity_StDev);
    cell_Affinity = IBM_Parameters.affinity*shifted_Distribution;
end

Define_Voxel_Domain;    %Define transitions in the voxel domain

NP_Evolution = zeros(IBM_Parameters.t_End/100,1);                                           %Number of associated nanoparticles over time
cell_Evolution = zeros(IBM_Parameters.t_End/100,1);                                         %Number of cells over time
NPs_per_Cell_Evolution = zeros(IBM_Parameters.t_End/100,IBM_Parameters.total_Cell_Spaces);  %Number of nanoparticles per cell over time
hist_Temp = cell(IBM_Parameters.t_End/100,1);                                               %Histogram of nanoparticles per cell
x_Range = cell(IBM_Parameters.t_End/100,1);                                                 %Range of histogram values from hist_Temp

Define_Distribution_Fits; %Define variables to store fit distributions

time_Save_Vec = 100:100:IBM_Parameters.t_End;                                               %Time points where data is stored
time_Save_Counter = 1;                                                                      %Counting variable of saved time points
time_Fit_Vec = 3600:3600:IBM_Parameters.t_End;                                              %Time points where distributions are fit
time_Fit_Counter = 1;                                                                       %Counting variable of saved fit points
sigma_Best_Vec = zeros(numel(time_Fit_Vec),1);                                              %Initialise best-fit heterogeneity vector
compound_Dist_Error = zeros(numel(time_Fit_Vec),1);                                         %Error in best-fit Poisson-lognormal distribution
compound_Dist_Error_Ideal = zeros(numel(time_Fit_Vec),1);                                   %Error in ideal Poisson-lognormal distribution
cell_Locations_Stored = cell(numel(time_Fit_Vec),1);                                        %Locations of cells at fit time points

if strcmpi(IBM_Parameters.cell_Cycle,'yes')
   number_G1_Cells = zeros(numel(time_Fit_Vec,1));          %Number of G1 phase cells  
   number_S_Cells = zeros(numel(time_Fit_Vec,1));           %Number of S phase cells
   number_G2_Cells = zeros(numel(time_Fit_Vec,1));          %Number of G2 phase cells
   number_M_Cells = zeros(numel(time_Fit_Vec,1));           %Number of M phase cells
end

if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
    mean_Affinity = zeros(numel(time_Fit_Vec),1);           %Effective mean affinity at fit time points
end

Define_Voxel_Indices;                                       %Define voxel indices for movement

Initialise_Cell_Cycle;

%% Voxel-based model simulation

while t < IBM_Parameters.t_End
    if sum(isnan(NP_Domain))~=0                             %Ensure no negative nanoparticles
        keyboard
    end
    
    number_Middle_NPs = sum(NP_Domain(middle_Indices));     %Number of nanoparticles in non-boundary voxels
    number_Top_NPs = sum(NP_Domain(top_Indices));           %Number of nanoparticles in top boundary voxels
    number_Bottom_NPs = sum(NP_Domain(bottom_Indices));     %Number of nanoparticles in bottom boundary voxels
    number_NPs = number_Middle_NPs+number_Top_NPs+number_Bottom_NPs; %Total number of nanoparticles
    
    propensity_Bind_Matrix(bottom_Indices) = (IBM_Parameters.D/IBM_Parameters.grid_Size^2+IBM_Parameters.V/(2*IBM_Parameters.grid_Size))*cell_Layer.*(cell_Affinity.*floor((cell_Capacity-NPs_per_Cell))./cell_Capacity); %Propensity for binding events with cell boundary

    Calculate_Event_Counts; %Calculate number of each event in time step
    
    Update_NP_Counts;       %Update the number of nanoparticles in each voxel
    
    %% Adaptive timestepping
    
    t = t + dt;             %Step time forward

    if min(min(NP_Domain)) < 0 || sum(NPs_per_Cell > cell_Layer.*cell_Capacity) > 0 %If number of nanoparticles below zero, or nanoparticles per cell above capacity, reduce time step and reset domain
        NP_Domain = NP_Domain_Temp; %Reset domain to before updated counts
        NPs_per_Cell = NPs_per_Cell - number_Bind_Events_Scaled(bottom_Indices); %Reset number of nanoparticles per cell
        t = t - dt;     %Reset time
        dt = dt/2;      %Temporarily halve timestep
    else
        if strcmpi(IBM_Parameters.cell_Cycle,'yes')
            Cell_Cycle_Update;  %Update cells via cell cycle
        end
        dt = IBM_Parameters.time_Step; %Reset timestep
    end
    
    %% Distribution fitting
    
    if t >= time_Save_Vec(time_Save_Counter) && max(NPs_per_Cell)>0 && t > 0 %Ensure number of nanoparticles per cell is positive and is a save time point
        counter = counter+1;
        
        x_Range{counter} = 0:max(NPs_per_Cell(cell_Locations)); %Calculate range of nanoparticles per cell
        
        NP_Evolution(counter) = mean(NPs_per_Cell);             %Store mean number of nanoparticles
        cell_Evolution(counter) = numel(cell_Locations);        %Store number of cells
        NPs_per_Cell_Evolution(counter,:) = NPs_per_Cell;       %Store number of nanoparticles per cell distribution
        hist_Temp{counter} = hist(NPs_per_Cell(cell_Locations),x_Range{counter})/numel(NPs_per_Cell(cell_Locations)); %Store histogram of number of nanoparticles per cell
        
        for i = 1:numel(IBM_Parameters.distribution_Fits)
            if  max(strcmpi(IBM_Parameters.distribution_Fits(i),["binomial";"poisson"]))
                evalc(sprintf('%s_Fit',IBM_Parameters.distribution_Fits(i)));               %Fit selected distributions
            end
        end
        
        if min(NPs_per_Cell(cell_Locations)) > 0
            
            for i = 1:numel(IBM_Parameters.distribution_Fits)
                if  max(strcmpi(IBM_Parameters.distribution_Fits(i),["lognormal";"normal";"negbin"])) 
                    evalc(sprintf('%s_Fit',IBM_Parameters.distribution_Fits(i)));           %Fit selected distributions
                end
            end
            
        end
        time_Save_Counter = time_Save_Counter + 1;
        t
    end
    
    if t >= time_Fit_Vec(time_Fit_Counter)                          
        hist_Data = hist(NPs_per_Cell(cell_Locations),0:max(NPs_per_Cell)); 
        sigma_Best_Vec(time_Fit_Counter) = lsqnonlin(@(sigma) compound_Distribution(NPs_per_Cell(cell_Locations),1000,0:max(NPs_per_Cell),sigma)-hist_Data,0.6,0,10,IBM_Parameters.opts); %Fit Poisson-lognormal distribution to the data
        tmpDiff = (compound_Distribution(NPs_per_Cell(cell_Locations),1000,0:max(NPs_per_Cell),sigma_Best_Vec(time_Fit_Counter))-hist_Data)/numel(NPs_per_Cell(cell_Locations));    %Error in best-fit Poisson-lognormal
        idealDiff = (compound_Distribution(NPs_per_Cell(cell_Locations),1000,0:max(NPs_per_Cell),IBM_Parameters.affinity_StDev)-hist_Data)/numel(NPs_per_Cell(cell_Locations));     %Error in ideal Poisson-lognormal
        compound_Dist_Error(time_Fit_Counter) = sqrt(sum(tmpDiff.^2));  %Error in best-fit Poisson-lognormal
        compound_Dist_Error_Ideal(time_Fit_Counter) = sqrt(sum(idealDiff.^2));  %Error in ideal Poisson-lognormal
        if strcmpi(IBM_Parameters.cell_Cycle,'yes')
            number_G1_Cells(time_Fit_Counter) = num_G1; %Store number of G1 cells
            number_S_Cells(time_Fit_Counter) = num_S;   %Store number of S cells
            number_G2_Cells(time_Fit_Counter) = num_G2; %Store number of G2 cells
            number_M_Cells(time_Fit_Counter) = num_M;   %Store number of M cells
        end
        if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
            mean_Affinity(time_Fit_Counter) = (num_G1*IBM_Parameters.G1_Affinity+num_S*IBM_Parameters.S_Affinity+num_G2*IBM_Parameters.G2_Affinity+num_M*IBM_Parameters.M_Affinity)/IBM_Parameters.number_of_Cells; %Store effective affinity
        end
        cell_Locations_Stored{time_Fit_Counter} = cell_Locations; %Store cell locations
        time_Fit_Counter = time_Fit_Counter+1;
    end
    
end

%% Calculate concentration profile in the z direction

vertical_Slice = zeros(IBM_Parameters.grid_Height,1);
for i = 1:IBM_Parameters.grid_Height
    vertical_Slice_Indices = i:IBM_Parameters.grid_Height:IBM_Parameters.grid_Height*IBM_Parameters.total_Cell_Spaces;
    vertical_Slice(i) = mean(NP_Domain(vertical_Slice_Indices));
end
