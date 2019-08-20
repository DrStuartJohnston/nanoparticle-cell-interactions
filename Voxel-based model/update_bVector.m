%% Update diagonal vector corresponding to boundary condition at cell layer for hybrid model

function bVector_End = update_bVector(PDE_Parameters,total_Associated,mean_dCell)

D = PDE_Parameters.rescaled_D;          %Diffusion coefficient
V = -PDE_Parameters.rescaled_V;         %Sedimentation coefficient
dt = PDE_Parameters.rescaled_dt;        %Time step
dx = PDE_Parameters.rescaled_dx;        %Node spacing                                                        
S = PDE_Parameters.rescaled_S;          %Confluency                                               

bVector_End = -(D*dt/dx^2+V*dt/(2*dx))*(-2*S*dx/D*mean_dCell-2*V*dx/D) + 1 + 2*D*dt/dx^2;