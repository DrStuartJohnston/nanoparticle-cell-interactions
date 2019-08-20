%% Define relevant inputs and calculate numerical solution for traditional dosage model

x = 0:PDE_Parameters.dx:PDE_Parameters.width;                                                           %Spatial domain

initial_Concentration_Value = 1.0;                                                                      %Value of initial condition
initial_Concentration_Width = 1.0*PDE_Parameters.width;                                                 %Width of initial condition

initial_Concentration = zeros(1,PDE_Parameters.num_Nodes);                                              %Initialise initial condition vector

for i = 1:PDE_Parameters.num_Nodes
    if x(i) <= initial_Concentration_Width
        initial_Concentration(i) = initial_Concentration_Value;                                         %Define initial condition
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
    PDE_Parameters.rescaled_x = x;                                             %Rescaled domain
    PDE_Parameters.rescaled_Final_Time = PDE_Parameters.final_Time;            %Rescaled end of experiment time
    PDE_Parameters.rescaled_dt = PDE_Parameters.dt;                            %Rescaled time step
    PDE_Parameters.rescaled_D = PDE_Parameters.D;                              %Rescaled diffusion coefficient
    PDE_Parameters.rescaled_V = PDE_Parameters.V;                              %Rescaled sedimentation coefficient
    PDE_Parameters.rescaled_dx = PDE_Parameters.dx;                            %Rescaled grid spacing
    PDE_Parameters.rescaled_S = PDE_Parameters.S;
    initial_Concentration = initial_Concentration*PDE_Parameters.particle_Concentration;
end

PDE_Parameters.initial = [initial_Concentration,0];                            %Pass in initial condition

time = 0:PDE_Parameters.rescaled_dt:PDE_Parameters.rescaled_Final_Time;        %Vector of solution time points

C = initial_Concentration';                                                    %Initial solution
total_Associated = zeros(numel(time),1);                                       %Vector of number of associated nanoparticles
[aVector,bVector,cVector] = define_ThomasAlg_Vectors(PDE_Parameters);          %Calculate vectors for Thomas' Algorithm 
opts = [];                                  

for tCounter = 1:numel(time)
    total_Associated(tCounter) = (sum(PDE_Parameters.initial)-sum(C))/sum(PDE_Parameters.initial)/PDE_Parameters.num_Cells*PDE_Parameters.initial_Particle_Count;   %Calculate number of associated nanoparticles
    bVector(end) = update_bVector_Single(PDE_Parameters,total_Associated(tCounter));    %Update vector for Thomas' Algorithm corresponding to boundary condition
    C = Thomas_Algorithm(aVector,bVector,cVector,C);    %Solve tridiagonal system via Thomas' algorithm           
end

final_Profile = C;                      %Final solution
scaled_Final_Profile = final_Profile;   %Final solution
uptake_Curve = total_Associated;        %Nanoparticle uptake curve

