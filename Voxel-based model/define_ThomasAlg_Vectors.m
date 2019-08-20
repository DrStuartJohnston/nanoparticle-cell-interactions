%% Function to define the tridiagonal vectors needed for Thomas' algorithm

function [aVector,bVector,cVector] = define_ThomasAlg_Vectors(PDE_Parameters)

aVector = zeros(PDE_Parameters.num_Nodes,1);    %Subdiagonal vector
bVector = zeros(PDE_Parameters.num_Nodes,1);    %Diagonal vector
cVector = zeros(PDE_Parameters.num_Nodes,1);    %Superdiagonal vector

D = PDE_Parameters.rescaled_D;                  %Diffusion coefficient
V = -PDE_Parameters.rescaled_V;                 %Sedimentation coefficient
dt = PDE_Parameters.rescaled_dt;                %Time step
dx = PDE_Parameters.rescaled_dx;                %Node spacing
        
bVector(2:end-1) = 1 + 2*D*dt/dx^2;             %Coefficients for P(x,t)
bVector(1) = 1 + 2*D*dt/dx^2 + (-D*dt/dx^2 + V*dt/dx)*(2*V*dx/D); %Coefficient for P(0,t)

aVector(2:end-1) = - D*dt/dx^2 + V*dt/(2*dx);   %Coefficients for P(x-h,t)
aVector(end) = -2*D*dt/dx^2;                    %Coefficient for P(L,t)
cVector(2:end-1) = - D*dt/dx^2 - V*dt/(2*dx);   %Coefficients for P(x+h,t)
cVector(1) = -2*D*dt/dx^2;                      %Coefficient for P(0,t)
