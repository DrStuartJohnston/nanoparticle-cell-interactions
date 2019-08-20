%% Initialise cell cycle parameters

total_Cycle_Time = IBM_Parameters.G1_Time+IBM_Parameters.S_Time+IBM_Parameters.G2_Time+IBM_Parameters.M_Time;   %Total time of cell cycle
G1_Phase_Proportion = IBM_Parameters.G1_Time/total_Cycle_Time;  %Proportion of cell cycle in G1 phase
S_Phase_Proportion = IBM_Parameters.S_Time/total_Cycle_Time;    %Proportion of cell cycle in S phase
G2_Phase_Proportion = IBM_Parameters.G2_Time/total_Cycle_Time;  %Proportion of cell cycle in G2 phase
M_Phase_Proportion = IBM_Parameters.M_Time/total_Cycle_Time;    %Proportion of cell cycle in M phase

%% Randomise phase of cell population
cell_Cycle_Rand = rand(IBM_Parameters.number_of_Cells,1);       
cell_Cycle_State = zeros(IBM_Parameters.number_of_Cells,1);

for i = 1:IBM_Parameters.number_of_Cells
    if cell_Cycle_Rand(i) < G1_Phase_Proportion
        cell_Cycle_State(i) = 1;    %1 corresponds to G1 phase
    elseif cell_Cycle_Rand(i) < S_Phase_Proportion+G1_Phase_Proportion
        cell_Cycle_State(i) = 2;    %2 corresponds to S phase
    elseif cell_Cycle_Rand(i) < G2_Phase_Proportion+S_Phase_Proportion+G1_Phase_Proportion
        cell_Cycle_State(i) = 3;    %3 corresponds to G2 phase
    else
        cell_Cycle_State(i) = 4;    %4 coresponds to M phase
    end
end

%% Define affinity and capacity parameters for each cell 
if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
    if strcmpi(IBM_Parameters.affinity_Variable,'yes')
        baseline_Affinity = lognrnd(log(1)-IBM_Parameters.capacity_StDev^2/2,IBM_Parameters.affinity_StDev,IBM_Parameters.total_Cell_Spaces,1); %Overall affinity distribution
        cell_Affinity = zeros(IBM_Parameters.total_Cell_Spaces,1);
        cell_Affinity(cell_Locations(cell_Cycle_State==1)) = IBM_Parameters.G1_Affinity*baseline_Affinity(cell_Locations(cell_Cycle_State==1)); %Affinity in G1 phase
        cell_Affinity(cell_Locations(cell_Cycle_State==2)) = IBM_Parameters.S_Affinity*baseline_Affinity(cell_Locations(cell_Cycle_State==2));  %Affinity in S phase
        cell_Affinity(cell_Locations(cell_Cycle_State==3)) = IBM_Parameters.G2_Affinity*baseline_Affinity(cell_Locations(cell_Cycle_State==3)); %Affinity in G2 phase
        cell_Affinity(cell_Locations(cell_Cycle_State==4)) = IBM_Parameters.M_Affinity*baseline_Affinity(cell_Locations(cell_Cycle_State==4));  %Affinity in M phase
        if strcmpi(IBM_Parameters.capacity_Variable,'yes') && strcmpi(IBM_Parameters.related_Capacity_Affinity,'yes')
            cell_Capacity = zeros(IBM_Parameters.total_Cell_Spaces,1);
            cell_Capacity = baseline_Affinity*IBM_Parameters.cell_Capacity; %Cell carrying capacity distribution
        end
    else
        cell_Affinity = zeros(IBM_Parameters.total_Cell_Spaces,1);
        cell_Affinity(cell_Locations(cell_Cycle_State==1)) = IBM_Parameters.G1_Affinity;    %Affinity in G1 phase
        cell_Affinity(cell_Locations(cell_Cycle_State==2)) = IBM_Parameters.S_Affinity;     %Affinity in S phase
        cell_Affinity(cell_Locations(cell_Cycle_State==3)) = IBM_Parameters.G2_Affinity;    %Affinity in G2 phase
        cell_Affinity(cell_Locations(cell_Cycle_State==4)) = IBM_Parameters.M_Affinity;     %Affinity in M phase
    end
end

num_G1 = sum(cell_Cycle_State==1);  %Number of G1 phase cells
num_S = sum(cell_Cycle_State==2);   %Number of S phase cells
num_G2 = sum(cell_Cycle_State==3);  %Number of G2 phase cells
num_M = sum(cell_Cycle_State==4);   %Number of M phase cells