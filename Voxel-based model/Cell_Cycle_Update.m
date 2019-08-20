G1_Locations = find(cell_Cycle_State==1);       %Locations of G1 cells
S_Locations = find(cell_Cycle_State==2);        %Locations of S cells
G2_Locations = find(cell_Cycle_State==3);       %Locations of G2 cells
M_Locations = find(cell_Cycle_State==4);        %Locations of M cells

num_G1 = numel(G1_Locations);                   %Number of G1 cells
num_S = numel(S_Locations);                     %Number of S cells
num_G2 = numel(G2_Locations);                   %Number of G2 cells
num_M = numel(M_Locations);                     %Number of M cells

propensity_G1 = 1/IBM_Parameters.G1_Time;       %Propensity of G1 transitions
propensity_S = 1/IBM_Parameters.S_Time;         %Propensity of S transitions
propensity_G2 = 1/IBM_Parameters.G2_Time;       %Propensity of G2 transitions
propensity_M = 1/IBM_Parameters.M_Time;         %Propensity of M transitions

events_G1 = poissrnd(propensity_G1*dt,num_G1,1);%Number of G1 transition events
events_S = poissrnd(propensity_S*dt,num_S,1);   %Number of S transition events
events_G2 = poissrnd(propensity_G2*dt,num_G2,1);%Number of G2 transition events
events_M = poissrnd(propensity_M*dt*(IBM_Parameters.total_Cell_Spaces-numel(cell_Locations))/IBM_Parameters.total_Cell_Spaces,num_M,1); %Number of M transition events


%% G1-S phase transition updates
for i = 1:num_G1
    if events_G1(i) > 0
        cell_Cycle_State(G1_Locations(i)) = 2;  %Define cell as S
        if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
            cell_Affinity(cell_Locations(G1_Locations(i))) = cell_Affinity(cell_Locations(G1_Locations(i)))/IBM_Parameters.G1_Affinity*IBM_Parameters.S_Affinity; %Update affinity
        end
    end
end

%% S-G2 phase transition updates
for i = 1:num_S
    if events_S(i) > 0
        cell_Cycle_State(S_Locations(i)) = 3; %Define cell as G2
        if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
            cell_Affinity(cell_Locations(S_Locations(i))) = cell_Affinity(cell_Locations(S_Locations(i)))/IBM_Parameters.S_Affinity*IBM_Parameters.G2_Affinity; %Update affinity
        end
    end
end

%% G2-M phase transition updates
for i = 1:num_G2
    if events_G2(i) > 0
        cell_Cycle_State(G2_Locations(i)) = 4; %Define cell as M
        if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
            cell_Affinity(cell_Locations(G2_Locations(i))) = cell_Affinity(cell_Locations(G2_Locations(i)))/IBM_Parameters.G2_Affinity*IBM_Parameters.M_Affinity; %Update affinity
        end
    end
end

%% M-G1 phase transition updates
for i = 1:num_M
    if events_M(i) > 0
        cell_Cycle_State(M_Locations(i)) = 1; %Define cell as G1
        potential_New_Site = find(cell_Layer==0); %Find new site for daughter cell
        rand_New_Site = ceil(numel(potential_New_Site)*rand);
        cell_Layer(potential_New_Site(rand_New_Site)) = 1;
        cell_Locations = [cell_Locations,potential_New_Site(rand_New_Site)];
        cell_Cycle_State = [cell_Cycle_State;1];
        NPs_per_Cell(potential_New_Site(rand_New_Site)) = floor(NPs_per_Cell(cell_Locations(M_Locations(i)))/2); %Split nanoparticle load
        if strcmpi(IBM_Parameters.capacity_Variable,'yes')
            cell_Capacity(potential_New_Site(rand_New_Site)) = cell_Capacity(cell_Locations(M_Locations(i))); %Inherit carrying capacity
        end
        if strcmpi(IBM_Parameters.affinity_Variable,'yes')
            cell_Affinity(potential_New_Site(rand_New_Site)) = cell_Affinity(cell_Locations(M_Locations(i))); %Inherit affinity
        end
        if strcmpi(IBM_Parameters.cell_Cycle_Dependent_Affinity,'yes')
            cell_Affinity(potential_New_Site(rand_New_Site)) = cell_Affinity(cell_Locations(M_Locations(i)))/IBM_Parameters.M_Affinity*IBM_Parameters.G1_Affinity; %Inherit affinity
            cell_Affinity(cell_Locations(M_Locations(i))) = cell_Affinity(cell_Locations(M_Locations(i)))/IBM_Parameters.M_Affinity*IBM_Parameters.G1_Affinity; %Inherit affinity
        end
        NPs_per_Cell(cell_Locations(M_Locations(i))) = ceil(NPs_per_Cell(cell_Locations(M_Locations(i)))/2); %Split nanoparticle load
    end
end

num_G1 = num_G1 - sum(events_G1>0) + 2*sum(events_M>0); %Number of G1 cells
num_S = num_S - sum(events_S>0) + sum(events_G1>0);     %Number of S cells
num_G2 = num_G2 - sum(events_G2>0) + sum(events_S>0);   %Number of G2 cells
num_M = num_M - sum(events_M>0) + sum(events_G2>0);     %Number of M cells

IBM_Parameters.number_of_Cells = IBM_Parameters.number_of_Cells + sum(events_M>0); %Total number of cells